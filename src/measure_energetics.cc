#include "measure_energetics.h"
#include "measure_energetics_style.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <utility>

// ---------------------------------------------------------------------- //
EnergeticsMeasure::EnergeticsMeasure(
    const Params& params,
    const Domain2D& domain,
    const Thermodynamics& thermodynamics,
    const FreeEnergy& free_energy,
    const TransportCoefficient& transport,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    thermodynamics_(thermodynamics),
    free_energy_(free_energy),
    transport_(transport),
    dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)),
    spectral_mask_(params, domain),
    num_order_parameters_(params.physics.num_order_parameters),
    num_fields_(params.physics.num_order_parameters + 3)
{
    if (num_order_parameters_ < 0) {
        throw std::runtime_error("EnergeticsMeasure requires a nonnegative number of order parameters.");
    }
    if (free_energy_.num_order_parameters() != num_order_parameters_) {
        throw std::runtime_error("EnergeticsMeasure received a free-energy model with inconsistent order-parameter count.");
    }

    auto cfg = std::dynamic_pointer_cast<const EnergeticsMeasureCommand>(command_);
    if (!cfg) {
        throw std::runtime_error("EnergeticsMeasure: invalid command type.");
    }

    nevery_ = cfg->nevery;
    file_ = cfg->file;

    k0_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    k2_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    k4_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);

    for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
        const std::size_t index = static_cast<std::size_t>(order_parameter);
        k0_coefficients_[index] = free_energy_.chemical_potential_k0_coefficient(order_parameter);
        k2_coefficients_[index] = free_energy_.chemical_potential_k2_coefficient(order_parameter);
        k4_coefficients_[index] = free_energy_.chemical_potential_k4_coefficient(order_parameter);
    }

    linear_pressure_coefficient_ = thermodynamics_.linear_pressure_coefficient();

}
// ---------------------------------------------------------------------- //
void EnergeticsMeasure::open_output_if_needed() {
    int open_ok = 1;

    if (domain_.rank() == 0 && !output_.is_open()) {
        output_.open(file_);
        if (!output_) {
            open_ok = 0;
        }
        else {
            output_ << "# step time total_energy kinetic_energy free_energy compressibility_energy\n";
            output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain_.comm());
    if (open_ok == 0) {
        throw std::runtime_error("EnergeticsMeasure: cannot open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::cell_area() const {
    return domain_.lx() * domain_.ly()
        / (static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global()));
}
// ---------------------------------------------------------------------- //
void EnergeticsMeasure::observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, int step, double time) {
    if (nevery_ <= 0 || step % nevery_ != 0) {
        return;
    }

    open_output_if_needed();

    const PhysicalStateBuffer* physical = nullptr;
    if (is_compressible_mode(dynamics_mode_) || thermodynamics_.has_physical_pressure() || free_energy_.has_physical_chemical_potential()) {
        physical = &workspace.ensure_physical(state, fft);
    }

    const double local_kinetic = compute_local_kinetic_energy(state, physical);
    const double local_free_energy = compute_local_free_energy(state, physical);
    const double local_compressibility = compute_local_compressibility_energy(state, physical);

    double local_values[3] = {local_kinetic, local_free_energy, local_compressibility};
    double global_values[3] = {0.0, 0.0, 0.0};

    MPI_Reduce(local_values, global_values, 3, MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

    int write_ok = 1;

    if (domain_.rank() == 0) {
        output_ << step
                << ' ' << time
                << ' ' << global_values[0] + global_values[1] + global_values[2]
                << ' ' << global_values[0]
                << ' ' << global_values[1]
                << ' ' << global_values[2]
                << '\n';

        if (!output_) {
            write_ok = 0;
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
    if (write_ok == 0) {
        throw std::runtime_error("EnergeticsMeasure: failed to write file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void EnergeticsMeasure::finalize() {
    if (output_.is_open()) {
        output_.close();
    }
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
    if (is_compressible_mode(dynamics_mode_)) {
        return compute_local_compressible_kinetic_energy(state, physical);
    }

    if (is_incompressible_mode(dynamics_mode_)) {
        return compute_local_incompressible_kinetic_energy(state, physical);
    }

    return 0.0;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
    double energy = 0.0;

    if (!is_compressible_mode(dynamics_mode_)) {
        return energy;
    }

    if (linear_pressure_coefficient_ != 0.0) {
        energy += compute_local_quadratic_compressibility_energy(state, physical);
    }

    if (thermodynamics_.has_physical_pressure()) {
        energy += compute_local_physical_compressibility_energy(state, physical);
    }

    return energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_free_energy(const State& state, const PhysicalStateBuffer* physical) const {
    double energy = 0.0;

    for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
        energy += compute_local_quadratic_free_energy(order_parameter, state, physical);
    }

    if (free_energy_.has_physical_chemical_potential()) {
        energy += compute_local_physical_free_energy(state, physical);
    }

    return energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::reference_density(const State& state) const {
    double local_rho0 = 0.0;
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    const Box2D& box = domain_.spectral_box();
    if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t index = static_cast<std::size_t>(0 - box.low[1]) * local_nkx + static_cast<std::size_t>(0 - box.low[0]);

        local_rho0 = state.rho_hat_data()[index].real() / grid_size;
    }

    double rho0 = 0.0;
    MPI_Allreduce(&local_rho0, &rho0, 1, MPI_DOUBLE, MPI_SUM, domain_.comm());

    return rho0;

}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_compressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
    (void)state;

    const Box2D& box = domain_.physical_box();
    const double dA = cell_area();

    double energy = 0.0;
    for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
        for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
            const double rho = physical->rho(gx, gy);
            const double jx = physical->jx(gx, gy);
            const double jy = physical->jy(gx, gy);
            energy += 0.5 * (jx * jx + jy * jy) / rho * dA;
        }
    }

    return energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_incompressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
    (void)physical;

    const double rho0 = reference_density(state);
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double prefactor = 0.5 * cell_area() / (grid_size * rho0);

    const Complex* jx = state.jx_hat_data();
    const Complex* jy = state.jy_hat_data();

    double energy = 0.0;
    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const std::size_t i = mode.index;
        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        energy += weight * (std::norm(jx[i]) + std::norm(jy[i]));
    }

    return prefactor * energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_quadratic_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
    (void)physical;
    double energy = 0.0;

    const double rho0 = reference_density(state);
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double prefactor = 0.5 * linear_pressure_coefficient_ * cell_area() / (grid_size * rho0);

    const Complex* rho = state.rho_hat_data();

    double spectral_sum = 0.0;
    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            continue;
        }

        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        spectral_sum += weight * std::norm(rho[mode.index]);
    }

    energy += prefactor * spectral_sum;

    return energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_physical_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
    (void)state;
    double energy = 0.0;

    const double dA = cell_area();
    const std::size_t local_size = domain_.physical_size();
    const double* rho = physical->rho_data();

    for (std::size_t i = 0; i < local_size; ++i) {
        energy += thermodynamics_.physical_thermodynamic_energy(rho[i]) * dA;
    }

    return energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_quadratic_free_energy(int order_parameter, const State& state, const PhysicalStateBuffer* physical) const {
    (void)physical;

    const std::size_t op = static_cast<std::size_t>(order_parameter);
    const double k0 = k0_coefficients_[op];
    const double k2_coeff = k2_coefficients_[op];
    const double k4_coeff = k4_coefficients_[op];

    if (k0 == 0.0 && k2_coeff == 0.0 && k4_coeff == 0.0) {
        return 0.0;
    }

    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double prefactor = 0.5 * cell_area() / grid_size;

    const Complex* psi = state.psi_hat_data(order_parameter);

    double energy = 0.0;
    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        const double k2 = mode.k2;
        const double coefficient = k0 + k2_coeff * k2 + k4_coeff * k2 * k2;
        const double weight = mode.gx == 0 ? 1.0 : 2.0;

        energy += weight * coefficient * std::norm(psi[mode.index]);
    }

    return prefactor * energy;
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::compute_local_physical_free_energy(const State& state, const PhysicalStateBuffer* physical) const {
    (void)state;

    double energy = 0.0;

    const double dA = cell_area();
    const std::size_t local_size = domain_.physical_size();

    std::vector<const double*> psi_physical(static_cast<std::size_t>(num_order_parameters_), nullptr);
    for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
        psi_physical[static_cast<std::size_t>(order_parameter)] = physical->psi_data(order_parameter);
    }

    std::vector<double> psi_point(static_cast<std::size_t>(num_order_parameters_), 0.0);

    for (std::size_t i = 0; i < local_size; ++i) {
        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            psi_point[static_cast<std::size_t>(order_parameter)] = psi_physical[static_cast<std::size_t>(order_parameter)][i];
        }

        energy += free_energy_.physical_free_energy_density(psi_point.data()) * dA;
    }

    return energy;
}
// ---------------------------------------------------------------------- //