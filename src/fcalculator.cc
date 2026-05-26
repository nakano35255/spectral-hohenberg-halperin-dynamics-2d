#include "fcalculator.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

// ---------------------------------------------------------------------- //
int FCalculator::temporary_field_capacity(const Params& params) {
    const int q = params.physics.num_order_parameters;

    if (q < 0) {
        throw std::runtime_error("FCalculator requires nonnegative num_order_parameters.");
    }

    return std::max(q, 4);
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
      num_order_parameters_(params.physics.num_order_parameters),
      local_spectral_size_(domain.spectral_size()),
      noise_rng_(),
      normal_noise_(0.0, 1.0),
      physical_rhs_plan_(make_physical_rhs_plan()),
      workspace_(domain, params, physical_rhs_plan_, temporary_field_capacity(params))
{
    if (num_order_parameters_ > 0 && free_energy_.has_physical_chemical_potential()) {
        const std::size_t rhs_size =
            static_cast<std::size_t>(num_order_parameters_) * local_spectral_size_;

        psi_nonlinear_rhs_.resize(rhs_size, Complex(0.0, 0.0));
        psi_point_.resize(static_cast<std::size_t>(num_order_parameters_), 0.0);
    }

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
bool FCalculator::is_incompressible_time_evolution() const {
    return params_.runtime.time_evolution_type == "euler/incompressible"
        || params_.runtime.time_evolution_type == "srk3/incompressible";
}
bool FCalculator::is_compressible_time_evolution() const {
    return params_.runtime.time_evolution_type == "euler/compressible"
        || params_.runtime.time_evolution_type == "srk3/compressible";
}
bool FCalculator::is_quiescent_time_evolution() const {
    return params_.runtime.time_evolution_type == "euler/quiescent"
        || params_.runtime.time_evolution_type == "srk3/quiescent";
}
// ---------------------------------------------------------------------- //
PhysicalRHSPlan FCalculator::make_physical_rhs_plan() const {
    PhysicalRHSPlan plan;
    plan.psi_det = make_psi_det_physical_request();
    plan.j_det = make_j_det_physical_request();
    return plan;
}
// ---------------------------------------------------------------------- //
PhysicalFieldRequestPlan FCalculator::make_psi_det_physical_request() const {
    PhysicalFieldRequestPlan request;

    if (num_order_parameters_ <= 0) {
        return request;
    }

    if (free_energy_.has_physical_chemical_potential()) {
        request.need_psi = true;
    }

    if (params_.fix.order_parameter_advection) {
        request.need_psi = true;
        request.need_j = true;
        request.need_velocity = true;

        if (is_compressible_time_evolution()) {
            request.need_rho = true;
        }
    }

    return request;
}
// ---------------------------------------------------------------------- //
PhysicalFieldRequestPlan FCalculator::make_j_det_physical_request() const {
    PhysicalFieldRequestPlan request;

    if (is_quiescent_time_evolution()) {
        return request;
    }

    const bool needs_physical_viscosity =
        transport_coefficient_.has_physical_shear_viscosity()
     || transport_coefficient_.has_physical_bulk_viscosity();

    if (needs_physical_viscosity) {
        throw std::runtime_error("physical viscosity is not implemented yet.");
    }

    if (is_compressible_time_evolution()) {
        request.need_rho = true;
        request.need_j = true;
        request.need_velocity = true;
    }

    if (params_.fix.momentum_advection) {
        request.need_j = true;
        request.need_velocity = true;
    }

    if (thermodynamics_.has_physical_pressure()) {
        if (!is_compressible_time_evolution()) {
            throw std::runtime_error("physical pressure requires compressible time evolution.");
        }
    }

    return request;
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
void FCalculator::ensure_psi_nonlinear_rhs(const State& current, double time) const {
    if (num_order_parameters_ == 0 || !free_energy_.has_physical_chemical_potential()) {
        return;
    }

    if (psi_nonlinear_rhs_ready_ && psi_nonlinear_state_ == &current && psi_nonlinear_time_ == time) {
        return;
    }

    psi_nonlinear_state_ = &current;
    psi_nonlinear_time_ = time;
    psi_nonlinear_rhs_ready_ = true;
    
    const std::size_t rhs_size = static_cast<std::size_t>(num_order_parameters_) * local_spectral_size_;
    std::fill(psi_nonlinear_rhs_.begin(), psi_nonlinear_rhs_.begin() + rhs_size, Complex(0.0, 0.0));

    const std::size_t local_physical_size = workspace_.local_physical_size();
    workspace_.ensure_physical_fields(current, time);
    const double* psi_physical = workspace_.physical_psi_data();
    double* mu_physical = workspace_.temporary_physical_fields(num_order_parameters_);
    Complex* mu_spectral = workspace_.temporary_spectral_fields(num_order_parameters_);

    for (std::size_t i = 0; i < local_physical_size; ++i) {
        for (int op = 0; op < num_order_parameters_; ++op) {
            psi_point_[static_cast<std::size_t>(op)] = psi_physical[static_cast<std::size_t>(op) * local_physical_size + i];
        }

        for (int op = 0; op < num_order_parameters_; ++op) {
            mu_physical[static_cast<std::size_t>(op) * local_physical_size + i] = free_energy_.physical_chemical_potential(op, psi_point_.data());
        }
    }

    workspace_.forward_many(num_order_parameters_, mu_physical, mu_spectral);

    const std::vector<double>& mobility = transport_coefficient_.order_parameter_mobility();

    for (int op = 0; op < num_order_parameters_; ++op) {
        const double mobility_op = mobility[static_cast<std::size_t>(op)];

        if (mobility_op == 0.0) {
            continue;
        }

        const std::size_t offset = static_cast<std::size_t>(op) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            psi_nonlinear_rhs_[offset + i] = -mobility_op * mode.k2 * mu_spectral[offset + i];
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::psi_det(int order_parameter, const State& current, Complex* out, double t) const {
    clear(out);

    const std::vector<double>& mobility = transport_coefficient_.order_parameter_mobility();
    const double mobility_op = mobility[static_cast<std::size_t>(order_parameter)];

    const bool order_parameter_advection = params_.fix.order_parameter_advection;

    if (mobility_op == 0.0 && !order_parameter_advection) {
        return;
    }

    const double mu_k0 = free_energy_.chemical_potential_k0_coefficient(order_parameter);
    const double mu_k2 = free_energy_.chemical_potential_k2_coefficient(order_parameter);
    const double mu_k4 = free_energy_.chemical_potential_k4_coefficient(order_parameter);

    const bool has_linear_mu = mu_k0 != 0.0 || mu_k2 != 0.0 || mu_k4 != 0.0;

    const Complex* psi = current.psi_hat_data(order_parameter);

    if (mobility_op != 0.0 && has_linear_mu) {
        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            const double k2 = mode.k2;
            const double k4 = k2 * k2;

            const double mu_coefficient = mu_k0 + mu_k2 * k2 + mu_k4 * k4;

            out[i] += -mobility_op * k2 * mu_coefficient * psi[i];
        }
    }

    if (mobility_op != 0.0 && free_energy_.has_physical_chemical_potential()) {
        ensure_psi_nonlinear_rhs(current, t);

        const std::size_t offset = static_cast<std::size_t>(order_parameter) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out[i] += psi_nonlinear_rhs_[offset + i];
        }
    }

    if (order_parameter_advection) {
        workspace_.ensure_physical_fields(current, t);

        const std::size_t local_physical_size = workspace_.local_physical_size();

        const double* psi_physical = workspace_.physical_psi_data() + static_cast<std::size_t>(order_parameter) * local_physical_size;

        const double* vx = workspace_.physical_vx_data();
        const double* vy = workspace_.physical_vy_data();

        double* flux_physical = workspace_.temporary_physical_fields(2);
        Complex* flux_spectral = workspace_.temporary_spectral_fields(2);

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            flux_physical[i] = psi_physical[i] * vx[i];
            flux_physical[local_physical_size + i] = psi_physical[i] * vy[i];
        }

        workspace_.forward_many(2, flux_physical, flux_spectral);

        const Complex* flux_x = flux_spectral;
        const Complex* flux_y = flux_spectral + local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out[i] += -Complex(0.0, 1.0) * (mode.kx * flux_x[i] + mode.ky * flux_y[i]);
        }
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::j_det(const State& current, Complex* out_jx, Complex* out_jy, double t) const {
    clear(out_jx);
    clear(out_jy);

    const double pressure_coefficient = thermodynamics_.linear_pressure_coefficient();
    const double eta = transport_coefficient_.shear_viscosity();
    const double zeta = transport_coefficient_.bulk_viscosity();

    const bool has_pressure = pressure_coefficient != 0.0;
    const bool has_viscosity = eta != 0.0 || zeta != 0.0;
    const bool momentum_advection = params_.fix.momentum_advection;

    if (!has_pressure && !has_viscosity && !momentum_advection) {
        return;
    }

    const Complex* rho = current.rho_hat_data();
    const Complex* jx = current.jx_hat_data();
    const Complex* jy = current.jy_hat_data();

    const bool compressible = is_compressible_time_evolution();

    Complex* velocity_spectral = nullptr;
    double rho0 = 1.0;

    if (has_viscosity) {
        if (compressible) {
            workspace_.ensure_physical_fields(current, t);
            velocity_spectral = workspace_.temporary_spectral_fields(2);

            workspace_.forward_many(2, workspace_.physical_vx_data(), velocity_spectral);
        }
        else {
            rho0 = workspace_.reference_density(current);
        }
    }

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double kx = mode.kx;
        const double ky = mode.ky;
        const double k2 = mode.k2;

        if (has_pressure) {
            out_jx[i] += -Complex(0.0, kx) * pressure_coefficient * rho[i];
            out_jy[i] += -Complex(0.0, ky) * pressure_coefficient * rho[i];
        }

        if (has_viscosity) {
            const Complex vx_hat = compressible ? velocity_spectral[i] : jx[i] / rho0;

            const Complex vy_hat = compressible ? velocity_spectral[local_spectral_size_ + i] : jy[i] / rho0;

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

    if (momentum_advection) {
        workspace_.ensure_physical_fields(current, t);

        const std::size_t local_physical_size = workspace_.local_physical_size();

        const double* jx_physical = workspace_.physical_jx_data();
        const double* jy_physical = workspace_.physical_jy_data();
        const double* vx = workspace_.physical_vx_data();
        const double* vy = workspace_.physical_vy_data();

        double* flux_physical = workspace_.temporary_physical_fields(4);
        Complex* flux_spectral = workspace_.temporary_spectral_fields(4);

        double* jx_vx = flux_physical;
        double* jx_vy = flux_physical + local_physical_size;
        double* jy_vx = flux_physical + 2 * local_physical_size;
        double* jy_vy = flux_physical + 3 * local_physical_size;

        for (std::size_t i = 0; i < local_physical_size; ++i) {
            jx_vx[i] = jx_physical[i] * vx[i];
            jx_vy[i] = jx_physical[i] * vy[i];
            jy_vx[i] = jy_physical[i] * vx[i];
            jy_vy[i] = jy_physical[i] * vy[i];
        }

        workspace_.forward_many(4, flux_physical, flux_spectral);

        const Complex* jx_vx_hat = flux_spectral;
        const Complex* jx_vy_hat = flux_spectral + local_spectral_size_;
        const Complex* jy_vx_hat = flux_spectral + 2 * local_spectral_size_;
        const Complex* jy_vy_hat = flux_spectral + 3 * local_spectral_size_;

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

            const std::size_t i = static_cast<std::size_t>(gy - box.low[1]) * local_nkx + local_kx0;
            const std::size_t ic = static_cast<std::size_t>(conjugate_gy - box.low[1]) * local_nkx + local_kx0;
            out[i] = value;
            out[ic] = std::conj(value);
        }
    }

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0) {
            continue;
        }

        if (mode.k2 == 0.0) {
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

            const std::size_t i = static_cast<std::size_t>(gy - box.low[1]) * local_nkx + local_kx0;
            const std::size_t ic = static_cast<std::size_t>(conjugate_gy - box.low[1]) * local_nkx + local_kx0;

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

