#include "measure_ave_profile.h"
#include "measure_ave_profile_style.h"

#include <algorithm>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
AveProfileMeasure::AveProfileMeasure(
     const Params& params,
     const Domain2D& domain,
     std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)),
    local_spectral_size_(domain.spectral_size()),
    local_physical_size_(domain.physical_size())
{
     auto cfg = std::dynamic_pointer_cast<const AveProfileMeasureCommand>(command_);
     if (!cfg) {
          throw std::runtime_error("AveProfileMeasure: invalid command type.");
     }

     axis_ = cfg->axis;
     nevery_ = cfg->nevery;
     nblock_ = cfg->nblock;
     file_ = cfg->file;
     average_mode_ = cfg->average_mode;
     targets_ = cfg->targets;

     profile_size_ = (axis_ == 0) ? domain_.nx_global() : domain_.ny_global();

     if (is_quiescent_mode(dynamics_mode_)) {
          for (const auto& target : targets_) {
               if (target.kind == AveProfileTargetKind::Vx || target.kind == AveProfileTargetKind::Vy) {
                    throw std::runtime_error("AveProfileMeasure vx/vy is not available in quiescent dynamics.");
               }
          }
     }

     spectral_slot_by_target_.assign(targets_.size(), -1);

     for (std::size_t target_index = 0; target_index < targets_.size(); ++target_index) {
          const AveProfileTarget& target = targets_[target_index];

          if (is_compressible_velocity_target(target)) {
               if (target.kind == AveProfileTargetKind::Vx && compressible_vx_slot_ < 0) {
                    compressible_vx_slot_ = num_compressible_velocity_fields_++;
               }
               if (target.kind == AveProfileTargetKind::Vy && compressible_vy_slot_ < 0) {
                    compressible_vy_slot_ = num_compressible_velocity_fields_++;
               }
               continue;
          }

          spectral_slot_by_target_[target_index] = num_spectral_output_fields_++;
     }

     const std::size_t total_spectral_size = static_cast<std::size_t>(num_spectral_output_fields_) * local_spectral_size_;

     spectral_block_sum_.assign(total_spectral_size, Complex(0.0, 0.0));
     spectral_running_sum_.assign(total_spectral_size, Complex(0.0, 0.0));

     if (num_compressible_velocity_fields_ > 0) {
          const std::size_t total_velocity_size = static_cast<std::size_t>(num_compressible_velocity_fields_) * local_physical_size_;

          compressible_velocity_block_sum_.assign(total_velocity_size, 0.0);
          compressible_velocity_running_sum_.assign(total_velocity_size, 0.0);
     }
}
// ---------------------------------------------------------------------- //
bool AveProfileMeasure::needs_incompressible_velocity() const {
     if (!is_incompressible_mode(dynamics_mode_)) {
          return false;
     }

     for (const auto& target : targets_) {
          if (target.kind == AveProfileTargetKind::Vx || target.kind == AveProfileTargetKind::Vy) {
               return true;
          }
     }
     return false;
}
// ---------------------------------------------------------------------- //
bool AveProfileMeasure::needs_order_parameter_flux() const {
     for (const auto& target : targets_) {
          if (target.kind == AveProfileTargetKind::OrderParameterFluxX ||
               target.kind == AveProfileTargetKind::OrderParameterFluxY) {
               return true;
          }
     }
     return false;
}
// ---------------------------------------------------------------------- //
bool AveProfileMeasure::needs_momentum_flux() const {
     for (const auto& target : targets_) {
          if (target.kind == AveProfileTargetKind::MomentumFluxXX ||
               target.kind == AveProfileTargetKind::MomentumFluxXY ||
               target.kind == AveProfileTargetKind::MomentumFluxYY) {
               return true;
          }
     }
     return false;
}
// ---------------------------------------------------------------------- //
bool AveProfileMeasure::is_compressible_velocity_target(const AveProfileTarget& target) const {
     if (!is_compressible_mode(dynamics_mode_)) {
          return false;
     }

     return target.kind == AveProfileTargetKind::Vx ||
            target.kind == AveProfileTargetKind::Vy;
}
// ---------------------------------------------------------------------- //
FluxRequest AveProfileMeasure::flux_request() const {
     FluxRequest request;
     request.order_parameter = needs_order_parameter_flux();
     request.momentum = needs_momentum_flux();
     return request;
}
// ---------------------------------------------------------------------- //
double AveProfileMeasure::reference_density(const State& state) {
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
          throw std::runtime_error("AveProfileMeasure encountered zero reference density.");
     }

     reference_density_ready_ = true;
     return reference_density_;
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::add_to_spectral_block(int slot, const Complex* values, double scale) {
     if (slot < 0) {
          throw std::runtime_error("AveProfileMeasure received an invalid spectral slot.");
     }

     Complex* dst = spectral_block_sum_.data() + static_cast<std::size_t>(slot) * local_spectral_size_;

     for (std::size_t i = 0; i < local_spectral_size_; ++i) {
          dst[i] += scale * values[i];
     }
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::add_to_compressible_velocity_block(int slot, const PhysicalStateBuffer& physical, int component) {
     if (slot < 0) {
          throw std::runtime_error("AveProfileMeasure received an invalid compressible velocity slot.");
     }

     double* dst = compressible_velocity_block_sum_.data() + static_cast<std::size_t>(slot) * local_physical_size_;

     const double* rho = physical.rho_data();
     const double* momentum = (component == 0) ? physical.jx_data() : physical.jy_data();

     for (std::size_t i = 0; i < local_physical_size_; ++i) {
          if (rho[i] == 0.0) {
               throw std::runtime_error("AveProfileMeasure compressible velocity encountered zero density.");
          }
          dst[i] += momentum[i] / rho[i];
     }
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::accumulate_sample(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux) {
     PhysicalStateBuffer* physical = nullptr;
     if (num_compressible_velocity_fields_ > 0) {
          physical = &workspace.ensure_physical(state, fft);

          if (compressible_vx_slot_ >= 0) {
               add_to_compressible_velocity_block(compressible_vx_slot_, *physical, 0);
          }
          if (compressible_vy_slot_ >= 0) {
               add_to_compressible_velocity_block(compressible_vy_slot_, *physical, 1);
          }
     }

     double rho0 = 0.0;
     if (needs_incompressible_velocity()) {
          rho0 = reference_density(state);
     }

     for (std::size_t target_index = 0; target_index < targets_.size(); ++target_index) {
          const int slot = spectral_slot_by_target_[target_index];
          if (slot < 0) {
               continue;
          }

          const AveProfileTarget& target = targets_[target_index];

          switch (target.kind) {
               case AveProfileTargetKind::Rho:
                    add_to_spectral_block(slot, state.rho_hat_data());
                    break;
               case AveProfileTargetKind::Psi:
                    add_to_spectral_block(slot, state.psi_hat_data(target.component));
                    break;
               case AveProfileTargetKind::Jx:
                    add_to_spectral_block(slot, state.jx_hat_data());
                    break;
               case AveProfileTargetKind::Jy:
                    add_to_spectral_block(slot, state.jy_hat_data());
                    break;
               case AveProfileTargetKind::Vx:
                    add_to_spectral_block(slot, state.jx_hat_data(), 1.0 / rho0);
                    break;
               case AveProfileTargetKind::Vy:
                    add_to_spectral_block(slot, state.jy_hat_data(), 1.0 / rho0);
                    break;
               case AveProfileTargetKind::OrderParameterFluxX:
                    add_to_spectral_block(slot, flux.order_parameter_flux_x_hat_data(target.component));
                    break;
               case AveProfileTargetKind::OrderParameterFluxY:
                    add_to_spectral_block(slot, flux.order_parameter_flux_y_hat_data(target.component));
                    break;
               case AveProfileTargetKind::MomentumFluxXX:
                    add_to_spectral_block(slot, flux.momentum_flux_xx_hat_data());
                    break;
               case AveProfileTargetKind::MomentumFluxXY:
                    add_to_spectral_block(slot, flux.momentum_flux_xy_hat_data());
                    break;
               case AveProfileTargetKind::MomentumFluxYY:
                    add_to_spectral_block(slot, flux.momentum_flux_yy_hat_data());
                    break;
          }
     }

     ++samples_in_block_;
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::finish_block(FourierTransform2D& fft, int step, double time) {
     if (samples_in_block_ <= 0) {
          throw std::runtime_error("AveProfileMeasure tried to finish an empty block.");
     }

     for (std::size_t i = 0; i < spectral_block_sum_.size(); ++i) {
          spectral_running_sum_[i] += spectral_block_sum_[i];
     }

     for (std::size_t i = 0; i < compressible_velocity_block_sum_.size(); ++i) {
          compressible_velocity_running_sum_[i] += compressible_velocity_block_sum_[i];
     }

     ++completed_blocks_;
     running_samples_ += samples_in_block_;

     const Complex* spectral_for_output = (average_mode_ == AveProfileAverageMode::Block) ? spectral_block_sum_.data() : spectral_running_sum_.data();
     const std::vector<double>* compressible_velocity_for_output = (average_mode_ == AveProfileAverageMode::Block) ? &compressible_velocity_block_sum_ : &compressible_velocity_running_sum_;
     const int output_samples = (average_mode_ == AveProfileAverageMode::Block) ? samples_in_block_ : running_samples_;

     const std::size_t total_profile_size = targets_.size() * static_cast<std::size_t>(profile_size_);
     std::vector<double> local_profiles(total_profile_size, 0.0);
     std::vector<double> global_profiles(total_profile_size, 0.0);
     const Box2D& box = domain_.physical_box();

     if (num_spectral_output_fields_ > 0) {
          const std::size_t total_physical_size = static_cast<std::size_t>(num_spectral_output_fields_) * local_physical_size_;

          std::vector<double> physical_for_output(total_physical_size, 0.0);

          fft.backward_many(num_spectral_output_fields_, spectral_for_output, physical_for_output.data());

          for (std::size_t target_index = 0; target_index < targets_.size(); ++target_index) {
               const int slot = spectral_slot_by_target_[target_index];
               if (slot < 0) {
                    continue;
               }

               const double* field = physical_for_output.data() + static_cast<std::size_t>(slot) * local_physical_size_;

               for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
                    for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                         const std::size_t local_index = domain_.physical_local_index(gx, gy);
                         const int bin = (axis_ == 0) ? gx : gy;
                         const std::size_t profile_index = target_index * static_cast<std::size_t>(profile_size_) + static_cast<std::size_t>(bin);

                         local_profiles[profile_index] += field[local_index];
                    }
               }
          }
     }

     if (num_compressible_velocity_fields_ > 0) {
          for (std::size_t target_index = 0; target_index < targets_.size(); ++target_index) {
               const AveProfileTarget& target = targets_[target_index];

               int slot = -1;
               if (target.kind == AveProfileTargetKind::Vx) {
                    slot = compressible_vx_slot_;
               } else if (target.kind == AveProfileTargetKind::Vy) {
                    slot = compressible_vy_slot_;
               } else {
                    continue;
               }

               if (slot < 0) {
                    continue;
               }

               const double* field = compressible_velocity_for_output->data() + static_cast<std::size_t>(slot) * local_physical_size_;

               for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
                    for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                         const std::size_t local_index = domain_.physical_local_index(gx, gy);
                         const int bin = (axis_ == 0) ? gx : gy;
                         const std::size_t profile_index = target_index * static_cast<std::size_t>(profile_size_) + static_cast<std::size_t>(bin);

                         local_profiles[profile_index] += field[local_index];
                    }
               }
          }
     }

     MPI_Reduce(local_profiles.data(), global_profiles.data(), static_cast<int>(global_profiles.size()), MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

     int write_ok = 1;
     if (domain_.rank() == 0) {
          const double transverse_grid_size = (axis_ == 0) ? static_cast<double>(domain_.ny_global()) : static_cast<double>(domain_.nx_global());
          const double normalizer = static_cast<double>(output_samples) * transverse_grid_size;

          std::vector<double> output_values(total_profile_size, 0.0);
          for (std::size_t i = 0; i < total_profile_size; ++i) {
               output_values[i] = global_profiles[i] / normalizer;
          }

          if (average_mode_ == AveProfileAverageMode::Block) {
               open_block_output_if_needed();
               write_rows(block_output_, step, time, output_samples, output_values);
               block_output_.flush();
               if (!block_output_) {
                    write_ok = 0;
               }
          } else {
               std::ofstream out(file_, std::ios::out | std::ios::trunc);
               if (!out) {
                    write_ok = 0;
               } else {
                    write_header(out);
                    out << std::scientific << std::setprecision(16);
                    write_rows(out, step, time, output_samples, output_values);
                    if (!out) {
                         write_ok = 0;
                    }
               }
          }
     }

     MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
     if (write_ok == 0) {
          throw std::runtime_error("AveProfileMeasure: failed to write file: " + file_);
     }

     std::fill(spectral_block_sum_.begin(), spectral_block_sum_.end(), Complex(0.0, 0.0));
     std::fill(compressible_velocity_block_sum_.begin(), compressible_velocity_block_sum_.end(), 0.0);
     samples_in_block_ = 0;
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::open_block_output_if_needed() {
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
          throw std::runtime_error("AveProfileMeasure: cannot open file: " + file_);
     }
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::write_header(std::ostream& out) const {
    out << "# measure ave/profile\n";
    out << "# axis " << (axis_ == 0 ? "x" : "y")
        << " nevery " << nevery_
        << " nblock " << nblock_
        << " average " << (average_mode_ == AveProfileAverageMode::Block ? "block" : "running")
        << "\n";
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::write_rows(std::ostream& out, int step, double time, int samples, const std::vector<double>& values) const {
     out << "# block " << completed_blocks_
          << " step " << step
          << " time " << time
          << " samples " << samples
          << "\n";

     out << "# columns coord";
     for (const auto& target : targets_) {
          out << ' ' << target.name;
     }
     out << "\n";

     const double spacing = (axis_ == 0) ? params_.grid.dx() : params_.grid.dy();

     for (int bin = 0; bin < profile_size_; ++bin) {
          const double coord = spacing * static_cast<double>(bin);
          out << coord;

          for (std::size_t target_index = 0; target_index < targets_.size(); ++target_index) {
               const std::size_t profile_index = target_index * static_cast<std::size_t>(profile_size_) + static_cast<std::size_t>(bin);
               out << ' ' << values[profile_index];
          }

          out << "\n";
     }

     out << "\n";
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) {
     ++block_step_;

     if (block_step_ % nevery_ == 0) {
          accumulate_sample(state, fft, workspace, flux);
     }

     if (block_step_ == nblock_) {
          finish_block(fft, step, time);
          block_step_ = 0;
     }
}
// ---------------------------------------------------------------------- //
void AveProfileMeasure::finalize() {
     if (block_output_.is_open()) {
          block_output_.close();
     }
}
// ---------------------------------------------------------------------- //
  
