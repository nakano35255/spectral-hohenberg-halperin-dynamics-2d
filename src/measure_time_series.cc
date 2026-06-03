#include "measure_time_series.h"
#include "measure_time_series_style.h"

#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <optional>
#include <vector>

// ---------------------------------------------------------------------- //
TimeSeriesMeasure::TimeSeriesMeasure(
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
    num_order_parameters_(params.physics.num_order_parameters)
{
    auto cfg = std::dynamic_pointer_cast<const TimeSeriesMeasureCommand>(command_);
    if (!cfg) {
        throw std::runtime_error("TimeSeriesMeasure: invalid command type.");
    }

    nevery_ = cfg->nevery;
    file_ = cfg->file;
    columns_ = cfg->columns;

    // for optional energy output
    (void)transport_;

    if (num_order_parameters_ < 0) {
        throw std::runtime_error("TimeSeriesMeasure requires a nonnegative number of order parameters.");
    }
    if (free_energy_.num_order_parameters() != num_order_parameters_) {
        throw std::runtime_error("TimeSeriesMeasure received a free-energy model with inconsistent order-parameter count.");
    }

    k0_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    k2_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    k4_coefficients_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    linear_pressure_coefficient_ = thermodynamics_.linear_pressure_coefficient();
    for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
        const std::size_t index = static_cast<std::size_t>(order_parameter);
        k0_coefficients_[index] = free_energy_.chemical_potential_k0_coefficient(order_parameter);
        k2_coefficients_[index] = free_energy_.chemical_potential_k2_coefficient(order_parameter);
        k4_coefficients_[index] = free_energy_.chemical_potential_k4_coefficient(order_parameter);
    }

}
// ---------------------------------------------------------------------- //
bool TimeSeriesMeasure::needs_order_parameter_flux() const {
    for (const TimeSeriesColumn& column : columns_) {
        if (column.kind == TimeSeriesColumnKind::OrderParameterFluxX ||
            column.kind == TimeSeriesColumnKind::OrderParameterFluxY) {
            return true;
        }
    }
    return false;
}
// ---------------------------------------------------------------------- //
bool TimeSeriesMeasure::needs_momentum_flux() const {
    for (const TimeSeriesColumn& column : columns_) {
        if (column.kind == TimeSeriesColumnKind::MomentumFluxXX ||
            column.kind == TimeSeriesColumnKind::MomentumFluxXY ||
            column.kind == TimeSeriesColumnKind::MomentumFluxYY) {
            return true;
        }
    }
    return false;
}
// ---------------------------------------------------------------------- //
bool TimeSeriesMeasure::needs_energy() const {
    for (const TimeSeriesColumn& column : columns_) {
        if (column.kind == TimeSeriesColumnKind::EnergyTotal ||
            column.kind == TimeSeriesColumnKind::EnergyKinetic ||
            column.kind == TimeSeriesColumnKind::EnergyFree ||
            column.kind == TimeSeriesColumnKind::EnergyCompressibility) {
            return true;
        }
    }
    return false;
}
// ---------------------------------------------------------------------- //
bool TimeSeriesMeasure::needs_spectral() const {
    for (const TimeSeriesColumn& column : columns_) {
        if (column.kind == TimeSeriesColumnKind::DensityL2 ||
            column.kind == TimeSeriesColumnKind::DensityGradientL2 ||
            column.kind == TimeSeriesColumnKind::MomentumL2 ||
            column.kind == TimeSeriesColumnKind::MomentumGradientL2 ||
            column.kind == TimeSeriesColumnKind::OrderParameterL2 ||
            column.kind == TimeSeriesColumnKind::OrderParameterGradientL2) {
            return true;
        }
    }
    return false;
}
// ---------------------------------------------------------------------- //
FluxRequest TimeSeriesMeasure::flux_request() const {
    FluxRequest request;
    request.order_parameter = needs_order_parameter_flux();
    request.momentum = needs_momentum_flux();
    return request;
}
// ---------------------------------------------------------------------- //
void TimeSeriesMeasure::open_output_if_needed() {
    int open_ok = 1;

    if (domain_.rank() == 0 && !output_.is_open()) {
        output_.open(file_);
        if (!output_) {
            open_ok = 0;
        } else {
            output_ << "# measure time_series\n";
            output_ << "# step time";
            for (const TimeSeriesColumn& column : columns_) {
                output_ << ' ' << column.name;
            }
            output_ << "\n";
            output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain_.comm());
    if (open_ok == 0) {
        throw std::runtime_error("TimeSeriesMeasure: cannot open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void TimeSeriesMeasure::observe(
    const State& state,
    FourierTransform2D& fft,
    MeasureWorkspace& workspace,
    const FluxBuffer& flux,
    int step,
    double time
) {
    if (nevery_ <= 0 || step % nevery_ != 0) {
        return;
    }

    open_output_if_needed();
    std::vector<double> local_values(columns_.size(), 0.0);

    // for optional energy output
    const bool output_energy = needs_energy();
    double local_kinetic = 0.0;
    double local_free_energy = 0.0;
    double local_compressibility = 0.0;
    if (output_energy) {
        const PhysicalStateBuffer* physical = nullptr;
        if (is_compressible_mode(dynamics_mode_) || thermodynamics_.has_physical_pressure() || free_energy_.has_physical_chemical_potential()) {
            physical = &workspace.ensure_physical(state, fft);
        }
        local_kinetic = compute_local_kinetic_energy(state, physical);
        local_free_energy = compute_local_free_energy(state, physical);
        local_compressibility = compute_local_compressibility_energy(state, physical);
    }

    for (std::size_t column_index = 0; column_index < columns_.size(); ++column_index) {
        const TimeSeriesColumn& column = columns_[column_index];

        switch (column.kind) {
            case TimeSeriesColumnKind::EnergyTotal:
                local_values[column_index] = local_kinetic + local_free_energy + local_compressibility;
                break;

            case TimeSeriesColumnKind::EnergyKinetic:
                local_values[column_index] = local_kinetic;
                break;

            case TimeSeriesColumnKind::EnergyFree:
                local_values[column_index] = local_free_energy;
                break;

            case TimeSeriesColumnKind::EnergyCompressibility:
                local_values[column_index] = local_compressibility;
                break;

            default:
                break;
        }
    }

    // for optional flux output
    const auto zero_index = zero_mode_index();
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    for (std::size_t column_index = 0; column_index < columns_.size(); ++column_index) {
        const TimeSeriesColumn& column = columns_[column_index];

        if (!zero_index) {
            continue;
        }

        switch (column.kind) {
            case TimeSeriesColumnKind::OrderParameterFluxX:
                local_values[column_index] =
                    flux.order_parameter_flux_x_hat_data(column.component)[*zero_index].real() / grid_size;
                break;

            case TimeSeriesColumnKind::OrderParameterFluxY:
                local_values[column_index] =
                    flux.order_parameter_flux_y_hat_data(column.component)[*zero_index].real() / grid_size;
                break;

            case TimeSeriesColumnKind::MomentumFluxXX:
                local_values[column_index] =
                    flux.momentum_flux_xx_hat_data()[*zero_index].real() / grid_size;
                break;

            case TimeSeriesColumnKind::MomentumFluxXY:
                local_values[column_index] =
                    flux.momentum_flux_xy_hat_data()[*zero_index].real() / grid_size;
                break;

            case TimeSeriesColumnKind::MomentumFluxYY:
                local_values[column_index] =
                    flux.momentum_flux_yy_hat_data()[*zero_index].real() / grid_size;
                break;

            default:
                break;
        }
    }

    // for optional L2 norm output
    for (std::size_t column_index = 0; column_index < columns_.size(); ++column_index) {
        const TimeSeriesColumn& column = columns_[column_index];

        switch (column.kind) {
            case TimeSeriesColumnKind::DensityL2:
                local_values[column_index] = compute_spectral_scalar_sum(state.rho_hat_data(), false);
                break;

            case TimeSeriesColumnKind::DensityGradientL2:
                local_values[column_index] = compute_spectral_scalar_sum(state.rho_hat_data(), true);
                break;

            case TimeSeriesColumnKind::MomentumL2:
                local_values[column_index] = compute_spectral_vector_sum(state.jx_hat_data(), state.jy_hat_data(), false);
                break;

            case TimeSeriesColumnKind::MomentumGradientL2:
                local_values[column_index] = compute_spectral_vector_sum(state.jx_hat_data(), state.jy_hat_data(), true);
                break;

            case TimeSeriesColumnKind::OrderParameterL2:
                local_values[column_index] = compute_spectral_scalar_sum(state.psi_hat_data(column.component), false);
                break;

            case TimeSeriesColumnKind::OrderParameterGradientL2:
                local_values[column_index] = compute_spectral_scalar_sum(state.psi_hat_data(column.component), true);
                break;

            default:
                break;
        }
    }

    std::vector<double> global_values(columns_.size(), 0.0);
    MPI_Reduce(local_values.data(), global_values.data(), static_cast<int>(local_values.size()), MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

    int write_ok = 1;
    if (domain_.rank() == 0) {
        output_ << step << ' ' << time;
        for (double value : global_values) {
            output_ << ' ' << value;
        }
        output_ << '\n';

        if (!output_) {
            write_ok = 0;
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
    if (write_ok == 0) {
        throw std::runtime_error("TimeSeriesMeasure: failed to write file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void TimeSeriesMeasure::finalize() {
    if (output_.is_open()) {
        output_.close();
    }
}
// ---------------------------------------------------------------------- //
double TimeSeriesMeasure::cell_area() const {
    return domain_.lx() * domain_.ly()
        / (static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global()));
}
// ---------------------------------------------------------------------- //
double TimeSeriesMeasure::compute_local_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
    if (is_compressible_mode(dynamics_mode_)) {
        return compute_local_compressible_kinetic_energy(state, physical);
    }

    if (is_incompressible_mode(dynamics_mode_)) {
        return compute_local_incompressible_kinetic_energy(state, physical);
    }

    return 0.0;
}
// ---------------------------------------------------------------------- //
double TimeSeriesMeasure::compute_local_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_free_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::reference_density(const State& state) const {
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
double TimeSeriesMeasure::compute_local_compressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_incompressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_quadratic_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_physical_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_quadratic_free_energy(int order_parameter, const State& state, const PhysicalStateBuffer* physical) const {
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
double TimeSeriesMeasure::compute_local_physical_free_energy(const State& state, const PhysicalStateBuffer* physical) const {
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
std::optional<std::size_t> TimeSeriesMeasure::zero_mode_index() const {
    const Box2D& local_box = domain_.spectral_box();

    if (local_box.low[0] <= 0 && 0 <= local_box.high[0] && local_box.low[1] <= 0 && 0 <= local_box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(local_box.size[0]);
        return static_cast<std::size_t>(0 - local_box.low[1]) * local_nkx
             + static_cast<std::size_t>(0 - local_box.low[0]);
    }

    return std::nullopt;
}
// ---------------------------------------------------------------------- //
double TimeSeriesMeasure::compute_spectral_scalar_sum(const Complex* field, bool with_k2) const {
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double prefactor = cell_area() / grid_size;

    double sum = 0.0;
    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            continue;
        }

        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        const double factor = with_k2 ? mode.k2 : 1.0;
        sum += weight * factor * std::norm(field[mode.index]);
    }

    return prefactor * sum;
}
// ---------------------------------------------------------------------- //
double TimeSeriesMeasure::compute_spectral_vector_sum(const Complex* field_x, const Complex* field_y, bool with_k2) const {
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double prefactor = cell_area() / grid_size;

    double sum = 0.0;
    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            continue;
        }

        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        const double factor = with_k2 ? mode.k2 : 1.0;
        const std::size_t i = mode.index;
        sum += weight * factor * (std::norm(field_x[i]) + std::norm(field_y[i]));
    }

    return prefactor * sum;
}
// ---------------------------------------------------------------------- //
