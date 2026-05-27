#include "fcalculator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
int FCalculator::temporary_field_capacity(const Params& params) {
    const int q = params.physics.num_order_parameters;

    if (q < 0) {
        throw std::runtime_error("FCalculator requires nonnegative num_order_parameters.");
    }

    return std::max(q, 7);
}
// ---------------------------------------------------------------------- //
FCalculator::FCalculator(
    const Params& params,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask,
    const Thermodynamics& thermodynamics,
    const FreeEnergy& free_energy,
    const TransportCoefficient& transport_coefficient
)
    : params_(params),
      domain_(domain),
      spectral_mask_(spectral_mask),
      thermodynamics_(thermodynamics),
      free_energy_(free_energy),
      transport_coefficient_(transport_coefficient),
      dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)),
      physical_field_requests_(PhysicalRHSFieldRequests::build(
          params,
          dynamics_mode_,
          free_energy,
          thermodynamics,
          transport_coefficient
      )),
      num_order_parameters_(params.physics.num_order_parameters),
      local_spectral_size_(domain.spectral_size()),
      workspace_(domain, params, dynamics_mode_, physical_field_requests_, temporary_field_capacity(params)),
      noise_rng_(),
      normal_noise_(0.0, 1.0)
{
    std::seed_seq noise_seed_sequence{
        static_cast<std::uint32_t>(params_.fix.noise.seed),
        static_cast<std::uint32_t>(domain_.rank())
    };
    noise_rng_.seed(noise_seed_sequence);
}
// ---------------------------------------------------------------------- //
void FCalculator::clear(Complex* out) const {
    std::fill(out, out + local_spectral_size_, Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
void FCalculator::rho_det(const State& current, Complex* out, double /*t*/) const {
    clear(out);

    const Complex* jx = current.jx_hat_data();
    const Complex* jy = current.jy_hat_data();

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        out[i] = -Complex(0.0, 1.0) * (mode.kx * jx[i] + mode.ky * jy[i]);
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::psi_det(int order_parameter, const State& current, Complex* out, double t) const {
    clear(out);

    const double mobility = transport_coefficient_.order_parameter_mobility()[static_cast<std::size_t>(order_parameter)];
    const bool has_advection = params_.fix.order_parameter_advection;
    const bool has_physical_mu = free_energy_.has_physical_chemical_potential();

    if (mobility == 0.0 && !has_advection) {
        return;
    }

    if (mobility != 0.0) {
        add_order_parameter_linear_term(order_parameter, mobility, current, out);
    }

    if (has_advection || has_physical_mu) {
        add_order_parameter_physical_terms(order_parameter, mobility, current, out, t);
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::j_det(const State& current, Complex* out_jx, Complex* out_jy, double t) const {
    clear(out_jx);
    clear(out_jy);

    if (is_quiescent_mode(dynamics_mode_)) {
        return;
    }

    const double pressure_coefficient = thermodynamics_.linear_pressure_coefficient();
    const double eta = transport_coefficient_.shear_viscosity();
    const double zeta = transport_coefficient_.bulk_viscosity();

    const bool has_linear_pressure = pressure_coefficient != 0.0;
    const bool has_physical_pressure = thermodynamics_.has_physical_pressure();
    const bool has_viscosity = eta != 0.0 || zeta != 0.0;
    const bool has_compressible_viscosity = has_viscosity && is_compressible_mode(dynamics_mode_);
    const bool has_incompressible_viscosity = has_viscosity && is_incompressible_mode(dynamics_mode_);
    const bool has_advection = params_.fix.momentum_advection;

    if (!has_linear_pressure && !has_physical_pressure && !has_viscosity && !has_advection) {
        return;
    }

    if (has_linear_pressure) {
        add_linear_pressure_term(pressure_coefficient, current, out_jx, out_jy);
    }

    if (has_incompressible_viscosity) {
        add_incompressible_viscous_term(eta, zeta, current, out_jx, out_jy);
    }

    if (has_compressible_viscosity) {
        add_compressible_viscous_term(eta, zeta, current, out_jx, out_jy, t);
    }

    if (has_physical_pressure || has_advection)
        add_momentum_physical_terms(eta, zeta, current, out_jx, out_jy, t
    );
}
// ---------------------------------------------------------------------- //
void FCalculator::add_order_parameter_linear_term(int order_parameter, double mobility, const State& current, Complex* out) const {
    const double mu_k0 = free_energy_.chemical_potential_k0_coefficient(order_parameter);
    const double mu_k2 = free_energy_.chemical_potential_k2_coefficient(order_parameter);
    const double mu_k4 = free_energy_.chemical_potential_k4_coefficient(order_parameter);

    if (mu_k0 == 0.0 && mu_k2 == 0.0 && mu_k4 == 0.0) {
        return;
    }

    const Complex* psi = current.psi_hat_data(order_parameter);

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double k2 = mode.k2;
        const double mu_coefficient = mu_k0 + mu_k2 * k2 + mu_k4 * k2 * k2;

        out[i] += -mobility * k2 * mu_coefficient * psi[i];
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::add_order_parameter_physical_terms(int order_parameter, double mobility, const State& current, Complex* out, double time) const {
    const bool has_physical_mu = mobility != 0.0 && free_energy_.has_physical_chemical_potential();
    const bool has_advection = params_.fix.order_parameter_advection;

    workspace_.ensure_physical_fields(current, time);

    const std::size_t local_physical_size = workspace_.local_physical_size();

    int nfields = 0;
    const int mu_field = has_physical_mu ? nfields++ : -1;
    const int flux_x_field = has_advection ? nfields++ : -1;
    const int flux_y_field = has_advection ? nfields++ : -1;

    double* physical_fields = workspace_.temporary_physical_fields(nfields);

    if (has_physical_mu) {
        const double* psi_physical = workspace_.physical_psi_data();
        double* mu_physical = physical_fields + static_cast<std::size_t>(mu_field) * local_physical_size;

        std::vector<double> psi_point(static_cast<std::size_t>(num_order_parameters_), 0.0);

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            for (int op = 0; op < num_order_parameters_; ++op) {
                psi_point[static_cast<std::size_t>(op)] = psi_physical[static_cast<std::size_t>(op) * local_physical_size + i];
            }
            mu_physical[i] = free_energy_.physical_chemical_potential(order_parameter, psi_point.data());
        }
    }

    if (has_advection) {
        const double* psi_physical = workspace_.physical_psi_data() + static_cast<std::size_t>(order_parameter) * local_physical_size;

        const double* vx = workspace_.physical_vx_data();
        const double* vy = workspace_.physical_vy_data();

        double* flux_x = physical_fields + static_cast<std::size_t>(flux_x_field) * local_physical_size;
        double* flux_y = physical_fields + static_cast<std::size_t>(flux_y_field) * local_physical_size;

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            flux_x[i] = psi_physical[i] * vx[i];
            flux_y[i] = psi_physical[i] * vy[i];
        }
    }

    Complex* spectral_fields = workspace_.temporary_spectral_fields(nfields);
    workspace_.forward_many(nfields, physical_fields, spectral_fields);

    if (has_physical_mu) {
        const Complex* mu_hat = spectral_fields + static_cast<std::size_t>(mu_field) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out[i] += -mobility * mode.k2 * mu_hat[i];
        }
    }

    if (has_advection) {
        const Complex* flux_x_hat = spectral_fields + static_cast<std::size_t>(flux_x_field) * local_spectral_size_;
        const Complex* flux_y_hat = spectral_fields + static_cast<std::size_t>(flux_y_field) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out[i] += -Complex(0.0, 1.0) * (mode.kx * flux_x_hat[i] + mode.ky * flux_y_hat[i]);
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::add_linear_pressure_term(double pressure_coefficient, const State& current, Complex* out_jx, Complex* out_jy) const {
    const Complex* rho = current.rho_hat_data();

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;

        out_jx[i] += -Complex(0.0, mode.kx) * pressure_coefficient * rho[i];
        out_jy[i] += -Complex(0.0, mode.ky) * pressure_coefficient * rho[i];
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::add_incompressible_viscous_term(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy) const {
    const double rho0 = workspace_.reference_density(current);

    const Complex* jx = current.jx_hat_data();
    const Complex* jy = current.jy_hat_data();

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double kx = mode.kx;
        const double ky = mode.ky;
        const double k2 = mode.k2;

        const Complex vx_hat = jx[i] / rho0;
        const Complex vy_hat = jy[i] / rho0;

        if (eta != 0.0) {
            out_jx[i] += -eta * k2 * vx_hat;
            out_jy[i] += -eta * k2 * vy_hat;
        }

        if (zeta != 0.0) {
            const Complex div_v = Complex(0.0, 1.0) * (kx * vx_hat + ky * vy_hat);
            out_jx[i] += zeta * Complex(0.0, kx) * div_v;
            out_jy[i] += zeta * Complex(0.0, ky) * div_v;
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::add_compressible_viscous_term(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy, double time) const {
    workspace_.ensure_physical_fields(current, time);

    const std::size_t local_physical_size = workspace_.local_physical_size();

    int nfields = 2;
    const int vx_field = 0;
    const int vy_field = 1;

    const double* vx = workspace_.physical_vx_data();
    const double* vy = workspace_.physical_vy_data();

    double* physical_fields = workspace_.temporary_physical_fields(nfields);
    double* vx_physical = physical_fields + static_cast<std::size_t>(vx_field) * local_physical_size;
    double* vy_physical = physical_fields + static_cast<std::size_t>(vy_field) * local_physical_size;

    std::copy(vx, vx + local_physical_size, vx_physical);
    std::copy(vy, vy + local_physical_size, vy_physical);

    Complex* spectral_fields = workspace_.temporary_spectral_fields(nfields);
    workspace_.forward_many(nfields, physical_fields, spectral_fields);

    const Complex* vx_hat = spectral_fields + static_cast<std::size_t>(vx_field) * local_spectral_size_;
    const Complex* vy_hat = spectral_fields + static_cast<std::size_t>(vy_field) * local_spectral_size_;

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double kx = mode.kx;
        const double ky = mode.ky;
        const double k2 = mode.k2;

        if (eta != 0.0) {
            out_jx[i] += -eta * k2 * vx_hat[i];
            out_jy[i] += -eta * k2 * vy_hat[i];
        }

        if (zeta != 0.0) {
            const Complex div_v = Complex(0.0, 1.0) * (kx * vx_hat[i] + ky * vy_hat[i]);
            out_jx[i] += zeta * Complex(0.0, kx) * div_v;
            out_jy[i] += zeta * Complex(0.0, ky) * div_v;
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::add_momentum_physical_terms(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy, double time) const {
    workspace_.ensure_physical_fields(current, time);

    const bool has_advection = params_.fix.momentum_advection;
    const bool has_physical_pressure = thermodynamics_.has_physical_pressure();

    const std::size_t local_physical_size = workspace_.local_physical_size();

    int nfields = 0;
    const int pressure_field = has_physical_pressure ? nfields++ : -1;
    const int jx_vx_field = has_advection ? nfields++ : -1;
    const int jx_vy_field = has_advection ? nfields++ : -1;
    const int jy_vx_field = has_advection ? nfields++ : -1;
    const int jy_vy_field = has_advection ? nfields++ : -1;

    double* physical_fields = workspace_.temporary_physical_fields(nfields);

    if (has_physical_pressure) {
        const double* rho = workspace_.physical_rho_data();
        double* pressure =
            physical_fields + static_cast<std::size_t>(pressure_field) * local_physical_size;

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            pressure[i] = thermodynamics_.physical_pressure(rho[i]);
        }
    }

    if (has_advection) {
        const double* jx = workspace_.physical_jx_data();
        const double* jy = workspace_.physical_jy_data();
        const double* vx = workspace_.physical_vx_data();
        const double* vy = workspace_.physical_vy_data();

        double* jx_vx = physical_fields + static_cast<std::size_t>(jx_vx_field) * local_physical_size;
        double* jx_vy = physical_fields + static_cast<std::size_t>(jx_vy_field) * local_physical_size;
        double* jy_vx = physical_fields + static_cast<std::size_t>(jy_vx_field) * local_physical_size;
        double* jy_vy = physical_fields + static_cast<std::size_t>(jy_vy_field) * local_physical_size;

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            jx_vx[i] = jx[i] * vx[i];
            jx_vy[i] = jx[i] * vy[i];
            jy_vx[i] = jy[i] * vx[i];
            jy_vy[i] = jy[i] * vy[i];
        }
    }

    Complex* spectral_fields = workspace_.temporary_spectral_fields(nfields);
    workspace_.forward_many(nfields, physical_fields, spectral_fields);

    if (has_physical_pressure) {
        const Complex* pressure_hat =
            spectral_fields + static_cast<std::size_t>(pressure_field) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out_jx[i] += -Complex(0.0, mode.kx) * pressure_hat[i];
            out_jy[i] += -Complex(0.0, mode.ky) * pressure_hat[i];
        }
    }

    if (has_advection) {
        const Complex* jx_vx_hat = spectral_fields + static_cast<std::size_t>(jx_vx_field) * local_spectral_size_;
        const Complex* jx_vy_hat = spectral_fields + static_cast<std::size_t>(jx_vy_field) * local_spectral_size_;
        const Complex* jy_vx_hat = spectral_fields + static_cast<std::size_t>(jy_vx_field) * local_spectral_size_;
        const Complex* jy_vy_hat = spectral_fields + static_cast<std::size_t>(jy_vy_field) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;

            out_jx[i] += -Complex(0.0, 1.0) * (mode.kx * jx_vx_hat[i] + mode.ky * jx_vy_hat[i]);

            out_jy[i] += -Complex(0.0, 1.0) * (mode.kx * jy_vx_hat[i] + mode.ky * jy_vy_hat[i]);
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::psi_sto(int order_parameter, const State& current, Complex* out) const {
    (void)current;

    clear(out);

    const std::vector<double>& mobility = transport_coefficient_.order_parameter_mobility();
    const double mobility_op = mobility[static_cast<std::size_t>(order_parameter)];
    const double kBT = params_.fix.noise.kBT;

    if (mobility_op == 0.0 || kBT == 0.0) {
        return;
    }

    const double active_grid_size = static_cast<double>(params_.grid.active_num_grid[0]) * static_cast<double>(params_.grid.active_num_grid[1]);
    const double volume = domain_.lx() * domain_.ly();
    const double prefactor = 2.0 * kBT * active_grid_size * active_grid_size / volume;

    const Box2D& box = domain_.spectral_box();
    const int ny = domain_.ny_global();

    if (box.low[0] == 0) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t local_kx0 = static_cast<std::size_t>(0 - box.low[0]);

        for (int gy = 0; gy < ny; ++gy) {
            if (!spectral_mask_.active(0, gy)) {
                continue;
            }

            const int conjugate_gy = (ny - gy) % ny;

            if (gy >= conjugate_gy) {
                continue;
            }

            const double ky = domain_.ky(gy);
            const double k2 = ky * ky;

            if (k2 == 0.0) {
                continue;
            }

            const double sigma = std::sqrt(0.5 * prefactor * mobility_op * k2);
            const Complex gaussian(normal_noise_(noise_rng_), normal_noise_(noise_rng_));
            const Complex value = sigma * gaussian;

            const std::size_t i =
                static_cast<std::size_t>(gy - box.low[1]) * local_nkx + local_kx0;
            const std::size_t ic =
                static_cast<std::size_t>(conjugate_gy - box.low[1]) * local_nkx + local_kx0;

            out[i] = value;
            out[ic] = std::conj(value);
        }
    }

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0 || mode.k2 == 0.0) {
            continue;
        }

        const double sigma = std::sqrt(0.5 * prefactor * mobility_op * mode.k2);
        const Complex gaussian(normal_noise_(noise_rng_), normal_noise_(noise_rng_));

        out[mode.index] = sigma * gaussian;
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::j_sto(const State& current, Complex* out_jx, Complex* out_jy) const {
    (void)current;

    clear(out_jx);
    clear(out_jy);

    const double eta = transport_coefficient_.shear_viscosity();
    const double zeta = transport_coefficient_.bulk_viscosity();
    const double kBT = params_.fix.noise.kBT;

    if (kBT == 0.0 || (eta == 0.0 && zeta == 0.0)) {
        return;
    }

    const double active_grid_size = static_cast<double>(params_.grid.active_num_grid[0]) * static_cast<double>(params_.grid.active_num_grid[1]);
    const double volume = domain_.lx() * domain_.ly();
    const double prefactor = 2.0 * kBT * active_grid_size * active_grid_size / volume;
    const double stress_sigma = std::sqrt(0.5 * prefactor);
    const double sqrt_eta = std::sqrt(eta);
    const double sqrt_zeta = std::sqrt(zeta);

    const Box2D& box = domain_.spectral_box();
    const int ny = domain_.ny_global();

    if (box.low[0] == 0) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t local_kx0 = static_cast<std::size_t>(0 - box.low[0]);

        for (int gy = 0; gy < ny; ++gy) {
            if (!spectral_mask_.active(0, gy)) {
                continue;
            }

            const int conjugate_gy = (ny - gy) % ny;

            if (gy >= conjugate_gy) {
                continue;
            }

            const double kx = 0.0;
            const double ky = domain_.ky(gy);
            const double k2 = ky * ky;

            if (k2 == 0.0) {
                continue;
            }

            const Complex xi1(normal_noise_(noise_rng_), normal_noise_(noise_rng_));
            const Complex xi2(normal_noise_(noise_rng_), normal_noise_(noise_rng_));
            const Complex xi3(normal_noise_(noise_rng_), normal_noise_(noise_rng_));

            const Complex sigma_xx = stress_sigma * (sqrt_zeta * xi1 + sqrt_eta * xi2);
            const Complex sigma_yy = stress_sigma * (sqrt_zeta * xi1 - sqrt_eta * xi2);
            const Complex sigma_xy = stress_sigma * sqrt_eta * xi3;

            const Complex value_x = Complex(0.0, 1.0) * (kx * sigma_xx + ky * sigma_xy);
            const Complex value_y = Complex(0.0, 1.0) * (kx * sigma_xy + ky * sigma_yy);

            const std::size_t i =
                static_cast<std::size_t>(gy - box.low[1]) * local_nkx + local_kx0;
            const std::size_t ic =
                static_cast<std::size_t>(conjugate_gy - box.low[1]) * local_nkx + local_kx0;

            out_jx[i] = value_x;
            out_jy[i] = value_y;
            out_jx[ic] = std::conj(value_x);
            out_jy[ic] = std::conj(value_y);
        }
    }

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0 || mode.k2 == 0.0) {
            continue;
        }

        const double kx = mode.kx;
        const double ky = mode.ky;

        const Complex xi1(normal_noise_(noise_rng_), normal_noise_(noise_rng_));
        const Complex xi2(normal_noise_(noise_rng_), normal_noise_(noise_rng_));
        const Complex xi3(normal_noise_(noise_rng_), normal_noise_(noise_rng_));

        const Complex sigma_xx = stress_sigma * (sqrt_zeta * xi1 + sqrt_eta * xi2);
        const Complex sigma_yy = stress_sigma * (sqrt_zeta * xi1 - sqrt_eta * xi2);
        const Complex sigma_xy = stress_sigma * sqrt_eta * xi3;

        out_jx[mode.index] = Complex(0.0, 1.0) * (kx * sigma_xx + ky * sigma_xy);
        out_jy[mode.index] = Complex(0.0, 1.0) * (kx * sigma_xy + ky * sigma_yy);
    }
}
// ---------------------------------------------------------------------- //