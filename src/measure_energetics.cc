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
    const Thermodynamics& thermodynamics,
    const FreeEnergy& free_energy,
    const TransportCoefficient& transport,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, std::move(command)),
    thermodynamics_(thermodynamics),
    free_energy_(free_energy),
    transport_(transport),
    dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)),
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
    psi_point_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);

    for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
        const std::size_t index = static_cast<std::size_t>(order_parameter);
        k0_coefficients_[index] = free_energy_.chemical_potential_k0_coefficient(order_parameter);
        k2_coefficients_[index] = free_energy_.chemical_potential_k2_coefficient(order_parameter);
        k4_coefficients_[index] = free_energy_.chemical_potential_k4_coefficient(order_parameter);
    }
}
// ---------------------------------------------------------------------- //
void EnergeticsMeasure::open_output_if_needed(const Domain2D& domain) {
    int open_ok = 1;

    if (domain.rank() == 0 && !output_.is_open()) {
        output_.open(file_);
        if (!output_) {
            open_ok = 0;
        }
        else {
            output_ << "# step time total_energy kinetic_energy free_energy compressibility_energy\n";
            output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain.comm());
    if (open_ok == 0) {
        throw std::runtime_error("EnergeticsMeasure: cannot open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
double EnergeticsMeasure::cell_area(const Domain2D& domain) const {
    return domain.lx() * domain.ly()
        / (static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global()));
}
// ---------------------------------------------------------------------- //
void EnergeticsMeasure::observe(
    const State& state,
    FourierTransform2D& fft,
    const Domain2D& domain,
    int step,
    double time
) {
    if (nevery_ <= 0 || step % nevery_ != 0) {
        return;
    }

    open_output_if_needed(domain);

    if (is_compressible_mode(dynamics_mode_) || thermodynamics_.linear_pressure_coefficient() || free_energy_.has_physical_chemical_potential()) {
        fft.backward_many(num_fields_, state.data(), physical.data());
    }

    local_kinetic = compute_local_kinetic_energy(state, domain);
    local_free_energy = compute_local_free_energy(state, domain);
    local_compressibility = compute_local_compressibility_energy(state, domain);

    double local_values[3] = {
        local_kinetic,
        local_free_energy,
        local_compressibility
    };
    double global_values[3] = {0.0, 0.0, 0.0};

    MPI_Reduce(local_values, global_values, 3, MPI_DOUBLE, MPI_SUM, 0, domain.comm());

    int write_ok = 1;

    if (domain.rank() == 0) {
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

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain.comm());
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
