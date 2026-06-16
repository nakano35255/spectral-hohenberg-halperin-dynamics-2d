#include "PASSIVE_SCALAR/measure_budget_spectrum.h"
#include "PASSIVE_SCALAR/measure_budget_spectrum_style.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------- //
namespace {
    constexpr double PI = 3.14159265358979323846;
}
// ---------------------------------------------------------------------- //
BudgetSpectrumMeasure::BudgetSpectrumMeasure(
    const Params& params,
    const Domain2D& domain,
    const FreeEnergy& free_energy,
    const TransportCoefficient& transport,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, std::move(command)),
    free_energy_(free_energy),
    transport_(transport),
    dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)),
    spectral_mask_(params, domain),
    num_order_parameters_(params.physics.num_order_parameters),
    local_physical_size_(domain.physical_size()),
    local_spectral_size_(domain.spectral_size())
{
    auto cfg = std::dynamic_pointer_cast<const BudgetSpectrumMeasureCommand>(command_);
    if (!cfg) {
        throw std::runtime_error("BudgetSpectrumMeasure: invalid command type.");
    }

    component_ = cfg->component;
    nevery_ = cfg->nevery;
    nblock_ = cfg->nblock;
    file_ = cfg->file;
    output_mode_ = cfg->output_mode;
    average_mode_ = cfg->average_mode;

    if (num_order_parameters_ < 0) {
        throw std::runtime_error("BudgetSpectrumMeasure requires a nonnegative number of order parameters.");
    }
    if (component_ < 0 || component_ >= num_order_parameters_) {
        throw std::runtime_error("BudgetSpectrumMeasure component index out of range.");
    }
    if (free_energy_.num_order_parameters() != num_order_parameters_) {
        throw std::runtime_error("BudgetSpectrumMeasure received an inconsistent free-energy model.");
    }
    if (transport_.num_order_parameters() != num_order_parameters_) {
        throw std::runtime_error("BudgetSpectrumMeasure received an inconsistent transport model.");
    }
    if (params_.fix.noise.order_parameter_enabled) {
        throw std::runtime_error("budget/spectrum does not include stochastic order-parameter noise.");
    }

    mobility_ = transport_.order_parameter_mobility()[static_cast<std::size_t>(component_)];
    k0_coefficient_ = free_energy_.chemical_potential_k0_coefficient(component_);
    k2_coefficient_ = free_energy_.chemical_potential_k2_coefficient(component_);
    k4_coefficient_ = free_energy_.chemical_potential_k4_coefficient(component_);
    has_physical_chemical_potential_ = free_energy_.has_physical_chemical_potential();

    has_transfer_ = params_.fix.order_parameter_advection;

    for (const SineForceFixConfig& force : params_.fix.sine_forces) {
        if (force.enabled && force.order_parameter_enabled && force.component == component_) {
            sine_forces_.push_back(force);
        }
    }

    for (const GradientForceFixConfig& force : params_.fix.gradient_forces) {
        if (force.enabled && force.component == component_) {
            gradient_forces_.push_back(force);
        }
    }

    has_production_ = !sine_forces_.empty() || !gradient_forces_.empty();

    const bool has_linear_chemical_potential = k0_coefficient_ != 0.0 || k2_coefficient_ != 0.0 || k4_coefficient_ != 0.0;

    has_dissipation_ = mobility_ != 0.0 && (has_linear_chemical_potential || has_physical_chemical_potential_);

    if ((has_transfer_ || !gradient_forces_.empty()) && is_quiescent_mode(dynamics_mode_)) {
        throw std::runtime_error("budget/spectrum cannot use velocity-dependent terms in quiescent dynamics.");
    }

    need_physical_velocity_ = has_transfer_ || (!gradient_forces_.empty() && is_compressible_mode(dynamics_mode_));

    need_velocity_hat_ = !gradient_forces_.empty();

    mode_block_sum_.assign(static_cast<std::size_t>(NUM_TERMS) * local_spectral_size_, 0.0);
    mode_running_sum_.assign(static_cast<std::size_t>(NUM_TERMS) * local_spectral_size_, 0.0);

    if (need_physical_velocity_) {
        velocity_physical_.assign(2 * local_physical_size_, 0.0);
    }

    if (need_velocity_hat_) {
        velocity_hat_.assign(2 * local_spectral_size_, Complex(0.0, 0.0));
    }

    if (has_transfer_) {
        transfer_flux_physical_.assign(2 * local_physical_size_, 0.0);
        transfer_flux_hat_.assign(2 * local_spectral_size_, Complex(0.0, 0.0));
    }

    if (has_dissipation_ && has_physical_chemical_potential_) {
        chemical_potential_physical_.assign(local_physical_size_, 0.0);
        chemical_potential_hat_.assign(local_spectral_size_, Complex(0.0, 0.0));
        psi_point_.assign(static_cast<std::size_t>(num_order_parameters_), 0.0);
    }

    dk_ = std::min(2.0 * PI / domain_.lx(), 2.0 * PI / domain_.ly());
    if (dk_ <= 0.0) {
        throw std::runtime_error("budget/spectrum requires positive domain lengths.");
    }

    initialize_shells();
}
// ---------------------------------------------------------------------- //
int BudgetSpectrumMeasure::shell_index(double k2) const {
    return static_cast<int>(std::floor(std::sqrt(k2) / dk_ + 0.5));
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::initialize_shells() {
    int local_max_shell = -1;

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.k2 == 0.0) {
            continue;
        }

        local_max_shell = std::max(local_max_shell, shell_index(mode.k2));
    }

    int global_max_shell = -1;
    MPI_Allreduce(&local_max_shell, &global_max_shell, 1, MPI_INT, MPI_MAX, domain_.comm());

    if (global_max_shell < 0) {
        throw std::runtime_error("budget/spectrum found no nonzero active spectral mode.");
    }

    nshells_ = global_max_shell + 1;
    shell_weight_counts_.assign(static_cast<std::size_t>(nshells_), 0.0);

    std::vector<double> local_counts(static_cast<std::size_t>(nshells_), 0.0);

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.k2 == 0.0) {
            continue;
        }

        const int shell = shell_index(mode.k2);
        const double weight = mode.gx == 0 ? 1.0 : 2.0;

        local_counts[static_cast<std::size_t>(shell)] += weight;
    }

    MPI_Allreduce(local_counts.data(), shell_weight_counts_.data(), nshells_, MPI_DOUBLE, MPI_SUM, domain_.comm());

    shell_block_sum_.assign(static_cast<std::size_t>(NUM_TERMS) * static_cast<std::size_t>(nshells_), 0.0);
    shell_running_sum_.assign(static_cast<std::size_t>(NUM_TERMS) * static_cast<std::size_t>(nshells_), 0.0);
}
// ---------------------------------------------------------------------- //
double BudgetSpectrumMeasure::reference_density(const State& state) {
    if (reference_density_ready_) {
        return reference_density_;
    }

    double local_rho0 = 0.0;
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    const Box2D& box = domain_.spectral_box();
    if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t index = static_cast<std::size_t>(0 - box.low[1]) * local_nkx +
                                    static_cast<std::size_t>(0 - box.low[0]);

        local_rho0 = state.rho_hat_data()[index].real() / grid_size;
    }

    MPI_Allreduce(&local_rho0, &reference_density_, 1, MPI_DOUBLE, MPI_SUM, domain_.comm());

    if (reference_density_ == 0.0) {
        throw std::runtime_error("budget/spectrum encountered zero reference density.");
    }

    reference_density_ready_ = true;
    return reference_density_;
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::prepare_velocity(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace) {
    if (!need_physical_velocity_ && !need_velocity_hat_) {
        return;
    }

    if (is_incompressible_mode(dynamics_mode_)) {
        const double rho0 = reference_density(state);

        if (need_physical_velocity_) {
            PhysicalStateBuffer& physical = workspace.ensure_physical(state, fft);

            const double* jx = physical.jx_data();
            const double* jy = physical.jy_data();

            for (std::size_t i = 0; i < local_physical_size_; ++i) {
                velocity_physical_[i] = jx[i] / rho0;
                velocity_physical_[local_physical_size_ + i] = jy[i] / rho0;
            }
        }

        if (need_velocity_hat_) {
            const Complex* jx_hat = state.jx_hat_data();
            const Complex* jy_hat = state.jy_hat_data();

            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                velocity_hat_[i] = jx_hat[i] / rho0;
                velocity_hat_[local_spectral_size_ + i] = jy_hat[i] / rho0;
            }
        }

        return;
    }

    if (is_compressible_mode(dynamics_mode_)) {
        PhysicalStateBuffer& physical = workspace.ensure_physical(state, fft);

        const double* rho = physical.rho_data();
        const double* jx = physical.jx_data();
        const double* jy = physical.jy_data();

        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            velocity_physical_[i] = jx[i] / rho[i];
            velocity_physical_[local_physical_size_ + i] = jy[i] / rho[i];
        }

        if (need_velocity_hat_) {
            fft.forward_many(2, velocity_physical_.data(), velocity_hat_.data());
        }

        return;
    }
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::prepare_transfer(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace) {
    if (!has_transfer_) {
        return;
    }

    PhysicalStateBuffer& physical = workspace.ensure_physical(state, fft);

    const double* psi = physical.psi_data(component_);
    const double* vx = velocity_physical_.data();
    const double* vy = velocity_physical_.data() + local_physical_size_;

    for (std::size_t i = 0; i < local_physical_size_; ++i) {
        transfer_flux_physical_[i] = psi[i] * vx[i];
        transfer_flux_physical_[local_physical_size_ + i] = psi[i] * vy[i];
    }

    fft.forward_many(2, transfer_flux_physical_.data(), transfer_flux_hat_.data());
}
// ---------------------------------------------------------------------- //
Complex BudgetSpectrumMeasure::transfer_term(const SpectralMode2D& mode) const {
    if (!has_transfer_) {
        return Complex(0.0, 0.0);
    }

    const Complex* flux_x_hat = transfer_flux_hat_.data();
    const Complex* flux_y_hat = transfer_flux_hat_.data() + local_spectral_size_;

    return -Complex(0.0, 1.0) * (mode.kx * flux_x_hat[mode.index] + mode.ky * flux_y_hat[mode.index]);
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::prepare_dissipation(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace) {
    if (!has_dissipation_ || !has_physical_chemical_potential_) {
        return;
    }

    PhysicalStateBuffer& physical = workspace.ensure_physical(state, fft);

    for (std::size_t i = 0; i < local_physical_size_; ++i) {
        for (int op = 0; op < num_order_parameters_; ++op) {
            psi_point_[static_cast<std::size_t>(op)] = physical.psi_data(op)[i];
        }

        chemical_potential_physical_[i] = free_energy_.physical_chemical_potential(component_, psi_point_.data());
    }

    fft.forward(chemical_potential_physical_.data(), chemical_potential_hat_.data());
}
// ---------------------------------------------------------------------- //
Complex BudgetSpectrumMeasure::dissipation_term(const SpectralMode2D& mode, const Complex& psi) const {
    if (!has_dissipation_) {
        return Complex(0.0, 0.0);
    }

    const double k2 = mode.k2;
    Complex value(0.0, 0.0);

    if (k0_coefficient_ != 0.0 || k2_coefficient_ != 0.0 || k4_coefficient_ != 0.0) {
        const double coefficient = k0_coefficient_ + k2_coefficient_ * k2 + k4_coefficient_ * k2 * k2;

        value += -mobility_ * k2 * coefficient * psi;
    }

    if (has_physical_chemical_potential_) {
        value += -mobility_ * k2 * chemical_potential_hat_[mode.index];
    }

    return value;
}
// ---------------------------------------------------------------------- //
Complex BudgetSpectrumMeasure::production_term(const SpectralMode2D& mode) const {
    if (!has_production_) {
        return Complex(0.0, 0.0);
    }

    Complex value(0.0, 0.0);

    const int ny = domain_.ny_global();
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    for (const SineForceFixConfig& force : sine_forces_) {
        const Complex positive_mode(0.0, -0.5 * force.amplitude * grid_size);
        const int gx_positive = (force.axis == 0) ? force.nk : 0;
        const int gy_positive = (force.axis == 0) ? 0 : force.nk;
        const int gy_negative = (ny - force.nk) % ny;

        if (mode.gx == gx_positive && mode.gy == gy_positive) {
            value += positive_mode;
        }

        if (force.axis == 1 && mode.gx == 0 && mode.gy == gy_negative) {
            value += std::conj(positive_mode);
        }
    }

    if (!gradient_forces_.empty()) {
        const Complex* vx_hat = velocity_hat_.data();
        const Complex* vy_hat = velocity_hat_.data() + local_spectral_size_;

        for (const GradientForceFixConfig& force : gradient_forces_) {
            const Complex* velocity = (force.direction == 0) ? vx_hat : vy_hat;
            value += -force.amplitude * velocity[mode.index];
        }
    }

    return value;
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::accumulate_sample(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace) {
    prepare_velocity(state, fft, workspace);
    prepare_transfer(state, fft, workspace);
    prepare_dissipation(state, fft, workspace);

    const Complex* psi_hat = state.psi_hat_data(component_);

    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double spectral_prefactor = 1.0 / (grid_size * grid_size);

    for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
        if (mode.k2 == 0.0) {
            continue;
        }

        const std::size_t i = mode.index;
        const Complex psi = psi_hat[i];

        const Complex transfer = transfer_term(mode);
        const Complex dissipation = dissipation_term(mode, psi);
        const Complex production = production_term(mode);

        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        const double factor = spectral_prefactor * weight;

        const std::array<double, NUM_TERMS> values = {
            factor * std::real(std::conj(psi) * transfer),
            factor * std::real(std::conj(psi) * dissipation),
            factor * std::real(std::conj(psi) * production)
        };

        const int shell = shell_index(mode.k2);

        for (int term = 0; term < NUM_TERMS; ++term) {
            mode_block_sum_[static_cast<std::size_t>(term) * local_spectral_size_ + i] += values[static_cast<std::size_t>(term)];
            shell_block_sum_[static_cast<std::size_t>(term) * static_cast<std::size_t>(nshells_) + static_cast<std::size_t>(shell)] += values[static_cast<std::size_t>(term)];
        }
    }

    ++samples_in_block_;
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::finish_block(int step, double time) {
    if (samples_in_block_ <= 0) {
        throw std::runtime_error("BudgetSpectrumMeasure tried to finish an empty block.");
    }

    for (std::size_t i = 0; i < mode_block_sum_.size(); ++i) {
        mode_running_sum_[i] += mode_block_sum_[i];
    }

    for (std::size_t i = 0; i < shell_block_sum_.size(); ++i) {
        shell_running_sum_[i] += shell_block_sum_[i];
    }

    ++completed_blocks_;
    running_samples_ += samples_in_block_;

    const int output_samples = (average_mode_ == BudgetSpectrumAverageMode::Block) ? samples_in_block_ : running_samples_;
    const std::vector<double>& mode_source = (average_mode_ == BudgetSpectrumAverageMode::Block) ? mode_block_sum_ : mode_running_sum_;
    const std::vector<double>& shell_source = (average_mode_ == BudgetSpectrumAverageMode::Block) ? shell_block_sum_ : shell_running_sum_;

    if (output_mode_ == BudgetSpectrumOutputMode::TwoD) {
        write_2d_output(step, time, output_samples, mode_source);
    } else {
        write_shell_output(step, time, output_samples, shell_source);
    }

    std::fill(mode_block_sum_.begin(), mode_block_sum_.end(), 0.0);
    std::fill(shell_block_sum_.begin(), shell_block_sum_.end(), 0.0);
    samples_in_block_ = 0;
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) {
    (void) flux;
    ++block_step_;

    if (block_step_ % nevery_ == 0) {
        accumulate_sample(state, fft, workspace);
    }

    if (block_step_ == nblock_) {
        finish_block(step, time);
        block_step_ = 0;
    }
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::finalize() {
    if (block_output_.is_open()) {
        block_output_.close();
    }
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::write_header(std::ostream& out) const {
    out << "# measure budget/spectrum\n";
    out << "# component " << component_
        << " nevery " << nevery_
        << " nblock " << nblock_
        << " mode " << (output_mode_ == BudgetSpectrumOutputMode::TwoD ? "2d" : "shell")
        << " average " << (average_mode_ == BudgetSpectrumAverageMode::Block ? "block" : "running")
        << " normalization per_volume"
        << "\n";
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::open_block_output_if_needed() {
    int open_ok = 1;

    if (domain_.rank() == 0 && !block_output_.is_open()) {
        block_output_.open(file_, std::ios::out);

        if (!block_output_) {
            open_ok = 0;
        } else {
            write_header(block_output_);
            block_output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain_.comm());

    if (open_ok == 0) {
        throw std::runtime_error("BudgetSpectrumMeasure: failed to open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::write_shell_output(int step, double time, int output_samples, const std::vector<double>& source) {
    if (average_mode_ == BudgetSpectrumAverageMode::Block) {
        open_block_output_if_needed();
    }

    std::vector<double> global_values;
    if (domain_.rank() == 0) {
        global_values.assign(source.size(), 0.0);
    }

    MPI_Reduce(source.data(), global_values.data(), static_cast<int>(source.size()), MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

    int write_ok = 1;

    if (domain_.rank() == 0) {
        std::ofstream running_output;
        std::ostream* out = nullptr;

        if (average_mode_ == BudgetSpectrumAverageMode::Block) {
            out = &block_output_;
        } else {
            running_output.open(file_, std::ios::out | std::ios::trunc);
            if (!running_output) {
                write_ok = 0;
            } else {
                write_header(running_output);
                running_output << std::scientific << std::setprecision(16);
                out = &running_output;
            }
        }

        if (write_ok == 1) {
            (*out) << "# block " << completed_blocks_
                   << " step " << step
                   << " time " << time
                   << " samples " << output_samples
                   << "\n";
            (*out) << "# k count transfer dissipation production total\n";

            for (int shell = 0; shell < nshells_; ++shell) {
                const double count = shell_weight_counts_[static_cast<std::size_t>(shell)];
                if (count == 0.0) {
                    continue;
                }

                const double transfer = global_values[static_cast<std::size_t>(TERM_TRANSFER) * nshells_ + shell] / static_cast<double>(output_samples);
                const double dissipation = global_values[static_cast<std::size_t>(TERM_DISSIPATION) * nshells_ + shell] / static_cast<double>(output_samples);
                const double production = global_values[static_cast<std::size_t>(TERM_PRODUCTION) * nshells_ + shell] / static_cast<double>(output_samples);
                const double total = transfer + dissipation + production;

                (*out) << static_cast<double>(shell) * dk_ << ' '
                       << count << ' '
                       << transfer << ' '
                       << dissipation << ' '
                       << production << ' '
                       << total << '\n';
            }

            (*out) << '\n';
            out->flush();

            if (!(*out)) {
                write_ok = 0;
            }
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());

    if (write_ok == 0) {
        throw std::runtime_error("BudgetSpectrumMeasure: failed to write file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void BudgetSpectrumMeasure::write_2d_output(int step, double time, int output_samples, const std::vector<double>& source) {
    if (average_mode_ == BudgetSpectrumAverageMode::Block) {
        open_block_output_if_needed();
    }

    const int local_count = static_cast<int>(source.size());

    std::vector<int> counts;
    if (domain_.rank() == 0) {
        counts.resize(static_cast<std::size_t>(domain_.size()));
    }

    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, domain_.comm());

    std::vector<int> displs;
    std::vector<double> gathered;

    if (domain_.rank() == 0) {
        displs.resize(static_cast<std::size_t>(domain_.size()), 0);

        int total_count = 0;
        for (int rank = 0; rank < domain_.size(); ++rank) {
            displs[static_cast<std::size_t>(rank)] = total_count;
            total_count += counts[static_cast<std::size_t>(rank)];
        }

        gathered.resize(static_cast<std::size_t>(total_count), 0.0);
    }

    MPI_Gatherv(source.data(), local_count, MPI_DOUBLE, gathered.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, domain_.comm());

    int write_ok = 1;

    if (domain_.rank() == 0) {
        const int nkx = domain_.nx_global() / 2 + 1;
        const int nky = domain_.ny_global();
        const std::size_t global_spectral_size = static_cast<std::size_t>(nkx) * static_cast<std::size_t>(nky);

        std::vector<double> global_values(static_cast<std::size_t>(NUM_TERMS) * global_spectral_size, 0.0);

        for (int rank = 0; rank < domain_.size(); ++rank) {
            const Box2D box = domain_.spectral_box_for_rank(rank);
            const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
            const std::size_t rank_spectral_size = static_cast<std::size_t>(box.size[0]) * static_cast<std::size_t>(box.size[1]);
            const std::size_t rank_offset = static_cast<std::size_t>(displs[static_cast<std::size_t>(rank)]);

            for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
                for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                    const std::size_t local_index = static_cast<std::size_t>(gy - box.low[1]) * local_nkx + static_cast<std::size_t>(gx - box.low[0]);
                    const std::size_t global_index = static_cast<std::size_t>(gy) * static_cast<std::size_t>(nkx) + static_cast<std::size_t>(gx);

                    for (int term = 0; term < NUM_TERMS; ++term) {
                        global_values[static_cast<std::size_t>(term) * global_spectral_size + global_index]
                         = gathered[rank_offset +static_cast<std::size_t>(term) * rank_spectral_size + local_index] / static_cast<double>(output_samples);
                    }
                }
            }
        }

        std::ofstream running_output;
        std::ostream* out = nullptr;

        if (average_mode_ == BudgetSpectrumAverageMode::Block) {
            out = &block_output_;
        } else {
            running_output.open(file_, std::ios::out | std::ios::trunc);
            if (!running_output) {
                write_ok = 0;
            } else {
                write_header(running_output);
                running_output << std::scientific << std::setprecision(16);
                out = &running_output;
            }
        }

        if (write_ok == 1) {
            (*out) << "# block " << completed_blocks_
                   << " step " << step
                   << " time " << time
                   << " samples " << output_samples
                   << "\n";
            (*out) << "# kx ky transfer dissipation production total\n";

            for (int gy = 0; gy < nky; ++gy) {
                for (int gx = 0; gx < nkx; ++gx) {
                    if (!spectral_mask_.active(gx, gy)) {
                        continue;
                    }

                    const double kx = domain_.kx(gx);
                    const double ky = domain_.ky(gy);
                    const double k2 = kx * kx + ky * ky;

                    if (k2 == 0.0) {
                        continue;
                    }

                    const std::size_t global_index = static_cast<std::size_t>(gy) * static_cast<std::size_t>(nkx) + static_cast<std::size_t>(gx);

                    const double transfer = global_values[static_cast<std::size_t>(TERM_TRANSFER) * global_spectral_size + global_index];
                    const double dissipation = global_values[static_cast<std::size_t>(TERM_DISSIPATION) * global_spectral_size + global_index];
                    const double production = global_values[static_cast<std::size_t>(TERM_PRODUCTION) * global_spectral_size + global_index];
                    const double total = transfer + dissipation + production;

                    (*out) << kx << ' '
                           << ky << ' '
                           << transfer << ' '
                           << dissipation << ' '
                           << production << ' '
                           << total << '\n';
                }
            }

            (*out) << '\n';
            out->flush();

            if (!(*out)) {
                write_ok = 0;
            }
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());

    if (write_ok == 0) {
        throw std::runtime_error("BudgetSpectrumMeasure: failed to write file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
