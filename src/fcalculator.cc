#include "fcalculator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double INV_TWO_POW_53 = 1.0 / 9007199254740992.0;

    std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
    }

    void hash_combine(std::uint64_t& seed, std::uint64_t value) {
        seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
    }

    double uniform_open_01(std::uint64_t key) {
        const std::uint64_t value = splitmix64(key);
        return (static_cast<double>(value >> 11) + 0.5) * INV_TWO_POW_53;
    }

    Complex standard_complex_gaussian(std::uint64_t key) {
        const double u1 = uniform_open_01(key);
        const double u2 = uniform_open_01(key + 0x9e3779b97f4a7c15ull);
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double angle = TWO_PI * u2;
        return Complex(radius * std::cos(angle), radius * std::sin(angle));
    }

    int signed_ky_index(int gy, int ny) {
        return (gy <= ny / 2) ? gy : gy - ny;
    }

    bool is_self_conjugate_x(int gx, int nx) {
        return gx == 0 || (nx % 2 == 0 && gx == nx / 2);
    }

    bool is_self_conjugate_y(int gy, int ny) {
        const int signed_ky = signed_ky_index(gy, ny);
        return signed_ky == 0 || (ny % 2 == 0 && std::abs(signed_ky) == ny / 2);
    }

    std::uint64_t noise_key(int seed, std::uint64_t call_count, int order_parameter, int gx, int gy) {
        std::uint64_t key = static_cast<std::uint64_t>(static_cast<std::int64_t>(seed));
        hash_combine(key, call_count);
        hash_combine(key, static_cast<std::uint64_t>(order_parameter));
        hash_combine(key, static_cast<std::uint64_t>(gx));
        hash_combine(key, static_cast<std::uint64_t>(gy));
        return key;
    }
}

// ---------------------------------------------------------------------- //
int FCalculator::temporary_field_capacity(const Params& params) {
    const int q = params.physics.num_order_parameters;

    if (q < 0) {
        throw std::runtime_error("FCalculator requires nonnegative num_order_parameters.");
    }

    // psi nonlinear: mu, psi vx, psi vy in future.
    // momentum nonlinear: pressure/advection buffers in future.
    return q;
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
      rhs_plan_(make_rhs_evaluation_plan()),
      workspace_(domain, params, temporary_field_capacity(params))
{
    if (num_order_parameters_ > 0 && free_energy_.has_physical_chemical_potential()) {
        const std::size_t rhs_size =
            static_cast<std::size_t>(num_order_parameters_) * local_spectral_size_;

        psi_nonlinear_rhs_.resize(rhs_size, Complex(0.0, 0.0));
        psi_point_.resize(static_cast<std::size_t>(num_order_parameters_), 0.0);
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::clear(Complex* out) const {
    std::fill(out, out + local_spectral_size_, Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
RHSEvaluationPlan FCalculator::make_rhs_evaluation_plan() const {
    RHSEvaluationPlan plan;
    plan.physical_fields = make_physical_field_plan();
    plan.dealias_rule = params_.grid.dealias_rule;
    return plan;
}
// ---------------------------------------------------------------------- //
PhysicalFieldPlan FCalculator::make_physical_field_plan() const {
    const bool needs_thermodynamics_physical = thermodynamics_.has_physical_pressure();
    const bool needs_free_energy_physical = free_energy_.has_physical_chemical_potential();
    const bool needs_transport_physical = transport_coefficient_.has_physical_shear_viscosity() || transport_coefficient_.has_physical_bulk_viscosity();

    if (!needs_thermodynamics_physical && !needs_free_energy_physical && !needs_transport_physical) {
        return PhysicalFieldPlan::None;
    }

    if (is_quiescent_time_evolution() && needs_free_energy_physical && !needs_thermodynamics_physical && !needs_transport_physical) {
        return PhysicalFieldPlan::PsiOnly;
    }

    return PhysicalFieldPlan::AllFields;
}
// ---------------------------------------------------------------------- //
bool FCalculator::is_quiescent_time_evolution() const {
    return params_.runtime.time_evolution_type == "euler/quiescent"
        || params_.runtime.time_evolution_type == "srk3/quiescent";
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
    const double* psi_physical = workspace_.physical_psi(current, time);
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
    if (mobility_op == 0.0) {
        return;
    }

    const double mu_k0 = free_energy_.chemical_potential_k0_coefficient(order_parameter);
    const double mu_k2 = free_energy_.chemical_potential_k2_coefficient(order_parameter);
    const double mu_k4 = free_energy_.chemical_potential_k4_coefficient(order_parameter);

    const bool has_linear_mu = mu_k0 != 0.0 || mu_k2 != 0.0 || mu_k4 != 0.0;

    const Complex* psi = current.psi_hat_data(order_parameter);

    if (has_linear_mu) {
        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            const double k2 = mode.k2;
            const double k4 = k2 * k2;

            const double mu_coefficient = mu_k0 + mu_k2 * k2 + mu_k4 * k4;

            out[i] += -mobility_op * k2 * mu_coefficient * psi[i];
        }
    }

    if (free_energy_.has_physical_chemical_potential()) {
        ensure_psi_nonlinear_rhs(current, t);

        const std::size_t offset = static_cast<std::size_t>(order_parameter) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t i = mode.index;
            out[i] += psi_nonlinear_rhs_[offset + i];
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

    const int compute_nx = domain_.nx_global();
    const int compute_ny = domain_.ny_global();
    const int active_nx = params_.grid.active_num_grid[0];
    const int active_ny = params_.grid.active_num_grid[1];
    const double active_grid_size = static_cast<double>(active_nx) * static_cast<double>(active_ny);
    const double volume = domain_.lx() * domain_.ly();
    const std::uint64_t call_count = psi_noise_call_count_++;

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;

        if (mode.k2 == 0.0) {
            out[i] = Complex(0.0, 0.0);
            continue;
        }

        const bool self_x = is_self_conjugate_x(mode.gx, compute_nx);
        const bool self_y = is_self_conjugate_y(mode.gy, compute_ny);
        const bool self_conjugate = self_x && self_y;

        int canonical_gy = mode.gy;
        bool use_conjugate = false;
        if (self_x && signed_ky_index(mode.gy, compute_ny) < 0) {
            canonical_gy = (compute_ny - mode.gy) % compute_ny;
            use_conjugate = true;
        }

        const std::uint64_t key = noise_key(
            params_.fix.noise.seed,
            call_count,
            order_parameter,
            mode.gx,
            canonical_gy
        );

        Complex gaussian = standard_complex_gaussian(key);
        if (use_conjugate) {
            gaussian = std::conj(gaussian);
        }
        if (self_conjugate) {
            gaussian = Complex(gaussian.real(), 0.0);
        }

        const double variance =
            2.0 * kBT * mobility_op * mode.k2 * active_grid_size * active_grid_size / volume;
        const double sigma = std::sqrt(self_conjugate ? variance : 0.5 * variance);
        out[i] = sigma * gaussian;
    }
}
// ---------------------------------------------------------------------- //
void FCalculator::j_det(const State& current, Complex* out_jx, Complex* out_jy, double /*t*/) const {
    clear(out_jx);
    clear(out_jy);

    const double pressure_coefficient = thermodynamics_.linear_pressure_coefficient();
    const double eta = transport_coefficient_.shear_viscosity();
    const double zeta = transport_coefficient_.bulk_viscosity();
    if (pressure_coefficient == 0.0 && eta == 0.0 && zeta == 0.0) {
        return;
    }

    const Complex* rho = current.rho_hat_data();
    const Complex* jx = current.jx_hat_data();
    const Complex* jy = current.jy_hat_data();

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double kx = mode.kx;
        const double ky = mode.ky;
        const double k2 = mode.k2;

        if (pressure_coefficient != 0.0) {
            out_jx[i] += -Complex(0.0, kx) * pressure_coefficient * rho[i];
            out_jy[i] += -Complex(0.0, ky) * pressure_coefficient * rho[i];
        }

        if (eta != 0.0) {
            out_jx[i] += -eta * k2 * jx[i];
            out_jy[i] += -eta * k2 * jy[i];
        }

        if (zeta != 0.0) {
            const Complex div_j = Complex(0.0, 1.0) * (kx * jx[i] + ky * jy[i]);
            out_jx[i] += zeta * Complex(0.0, kx) * div_j;
            out_jy[i] += zeta * Complex(0.0, ky) * div_j;
        }
    }
}
// ---------------------------------------------------------------------- //
