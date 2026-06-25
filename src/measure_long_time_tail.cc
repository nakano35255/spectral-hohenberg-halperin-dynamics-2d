#include "measure_long_time_tail.h"
#include "measure_long_time_tail_style.h"

#include <algorithm>
#include <complex>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
LongTimeTailMeasure::LongTimeTailMeasure(
     const Params& params,
     const Domain2D& domain,
     std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    spectral_mask_(params, domain),
    local_spectral_size_(domain.spectral_size())
{
     auto cfg = std::dynamic_pointer_cast<const LongTimeTailMeasureCommand>(command_);
     if (!cfg) {
          throw std::runtime_error("LongTimeTailMeasure: invalid command type.");
     }

     nevery_ = cfg->nevery;
     nblock_ = cfg->nblock;
     file_ = cfg->file;
     average_mode_ = cfg->average_mode;
     targets_ = cfg->targets;
     cross_ = cfg->cross;

     nfields_ = static_cast<int>(targets_.size());
     nlags_ = nblock_ / nevery_;
     npairs_ = cross_ ? nfields_ * (nfields_ + 1) / 2 : nfields_;

     block_fields_.assign(
          static_cast<std::size_t>(nlags_) *
          static_cast<std::size_t>(nfields_) *
          local_spectral_size_,
          Complex(0.0, 0.0)
     );

     block_sum_.assign(
          static_cast<std::size_t>(npairs_) * static_cast<std::size_t>(nlags_),
          0.0
     );
     running_sum_.assign(block_sum_.size(), 0.0);
}
// ---------------------------------------------------------------------- //
std::size_t LongTimeTailMeasure::sample_field_offset(int sample_index, int field_index) const {
     return (
          static_cast<std::size_t>(sample_index) * static_cast<std::size_t>(nfields_) +
          static_cast<std::size_t>(field_index)
     ) * local_spectral_size_;
}
// ---------------------------------------------------------------------- //
std::size_t LongTimeTailMeasure::pair_lag_index(int pair_index, int lag_index) const {
     return (
          static_cast<std::size_t>(pair_index) * static_cast<std::size_t>(nlags_) +
          static_cast<std::size_t>(lag_index)
     );
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::store_current_sample(const State& state, int sample_index) {
     for (int field_index = 0; field_index < nfields_; ++field_index) {
          const LongTimeTailTarget& target = targets_[static_cast<std::size_t>(field_index)];
          Complex* dst = block_fields_.data() + sample_field_offset(sample_index, field_index);
          const Complex* src = nullptr;

          switch (target.kind) {
               case LongTimeTailTargetKind::Rho:
                    src = state.rho_hat_data();
                    break;
               case LongTimeTailTargetKind::Jx:
                    src = state.jx_hat_data();
                    break;
               case LongTimeTailTargetKind::Jy:
                    src = state.jy_hat_data();
                    break;
               case LongTimeTailTargetKind::Psi:
                    src = state.psi_hat_data(target.component);
                    break;
          }

          std::copy(src, src + local_spectral_size_, dst);
     }
}
// ---------------------------------------------------------------------- //
double LongTimeTailMeasure::compute_local_correlation(const Complex* later, const Complex* earlier) const {
     double sum = 0.0;

     for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
          if (mode.gx == 0 && mode.gy == 0) {
               continue;
          }

          const double r2c_weight = (mode.gx == 0) ? 1.0 : 2.0;
          sum += r2c_weight * (later[mode.index] * std::conj(earlier[mode.index])).real();
     }

     const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

     return sum / (grid_size * grid_size);
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::accumulate_auto_correlations(std::vector<double>& local_sum) const {
     for (int lag_index = 0; lag_index < nlags_; ++lag_index) {
          const int count = nlags_ - lag_index;

          for (int origin = 0; origin < count; ++origin) {
               for (int field_index = 0; field_index < nfields_; ++field_index) {
                    const Complex* later = block_fields_.data() + sample_field_offset(origin + lag_index, field_index);
                    const Complex* earlier = block_fields_.data() + sample_field_offset(origin, field_index);

                    const int pair_index = field_index;
                    local_sum[pair_lag_index(pair_index, lag_index)] += compute_local_correlation(later, earlier);
               }
          }
     }
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::accumulate_cross_correlations(std::vector<double>& local_sum) const {
     for (int lag_index = 0; lag_index < nlags_; ++lag_index) {
          const int count = nlags_ - lag_index;
          int pair_index = 0;

          for (int later_field = 0; later_field < nfields_; ++later_field) {
               for (int earlier_field = later_field; earlier_field < nfields_; ++earlier_field) {
                    for (int origin = 0; origin < count; ++origin) {
                         const Complex* later = block_fields_.data() + sample_field_offset(origin + lag_index, later_field);
                         const Complex* earlier = block_fields_.data() + sample_field_offset(origin, earlier_field);

                         local_sum[pair_lag_index(pair_index, lag_index)] += compute_local_correlation(later, earlier);
                    }

                    ++pair_index;
               }
          }
     }
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::finish_block(int step, double time) {
     std::vector<double> local_sum(block_sum_.size(), 0.0);

     if (cross_) {
          accumulate_cross_correlations(local_sum);
     } else {
          accumulate_auto_correlations(local_sum);
     }

     MPI_Allreduce(local_sum.data(), block_sum_.data(), static_cast<int>(block_sum_.size()), MPI_DOUBLE, MPI_SUM, domain_.comm());

     for (std::size_t i = 0; i < running_sum_.size(); ++i) {
          running_sum_[i] += block_sum_[i];
     }

     ++completed_blocks_;

     const std::vector<double>& sums = (average_mode_ == LongTimeTailAverageMode::Block) ? block_sum_ : running_sum_;

     if (average_mode_ == LongTimeTailAverageMode::Block) {
          open_block_output_if_needed();
     }

     int write_ok = 1;
     if (domain_.rank() == 0) {
          if (average_mode_ == LongTimeTailAverageMode::Block) {
               write_rows(block_output_, step, time, sums);
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
                    write_rows(out, step, time, sums);
                    if (!out) {
                         write_ok = 0;
                    }
               }
          }
     }

     MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
     if (write_ok == 0) {
          throw std::runtime_error("LongTimeTailMeasure: failed to write file: " + file_);
     }

     reset_block();
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::reset_block() {
     std::fill(block_sum_.begin(), block_sum_.end(), 0.0);
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::open_block_output_if_needed() {
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
          throw std::runtime_error("LongTimeTailMeasure: cannot open file: " + file_);
     }
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::write_header(std::ostream& out) const {
     out << "# measure long_time_tail\n";
     out << "# nevery " << nevery_
         << " nblock " << nblock_
         << " average " << (average_mode_ == LongTimeTailAverageMode::Block ? "block" : "running")
         << " cross " << (cross_ ? "on" : "off")
         << "\n";

     out << "# columns nsamples tau";

     if (cross_) {
          for (int later_field = 0; later_field < nfields_; ++later_field) {
               for (int earlier_field = later_field; earlier_field < nfields_; ++earlier_field) {
                    out << ' ' << targets_[static_cast<std::size_t>(later_field)].name
                        << targets_[static_cast<std::size_t>(earlier_field)].name;
               }
          }
     } else {
          for (int field_index = 0; field_index < nfields_; ++field_index) {
               out << ' ' << targets_[static_cast<std::size_t>(field_index)].name
                   << targets_[static_cast<std::size_t>(field_index)].name;
          }
     }

     out << "\n";
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::write_rows(std::ostream& out, int, double, const std::vector<double>& sums) const {
     if (average_mode_ == LongTimeTailAverageMode::Block) {
          out << "# block " << completed_blocks_ << "\n";
     }

     for (int lag_index = 0; lag_index < nlags_; ++lag_index) {
          const double tau = static_cast<double>(lag_index * nevery_) * params_.runtime.dt;

          const long long block_count = static_cast<long long>(nlags_ - lag_index);
          const long long count = (average_mode_ == LongTimeTailAverageMode::Block) ? block_count : static_cast<long long>(completed_blocks_) * block_count;

          out << count << ' ' << tau;

          for (int pair_index = 0; pair_index < npairs_; ++pair_index) {
               const double value = sums[pair_lag_index(pair_index, lag_index)] / static_cast<double>(count);

               out << ' ' << value;
          }

          out << '\n';
     }

     out << '\n';
}
// ---------------------------------------------------------------------- //
void LongTimeTailMeasure::observe(const State& state, FourierTransform2D&, MeasureWorkspace&, const FluxBuffer&, int step, double time) {
     ++block_step_;

     if (block_step_ % nevery_ == 0) {
          const int sample_index = block_step_ / nevery_ - 1;
          store_current_sample(state, sample_index);
     }

     if (block_step_ == nblock_) {
          finish_block(step, time);
          block_step_ = 0;
     }
}

void LongTimeTailMeasure::finalize() {
     if (block_output_.is_open()) {
          block_output_.close();
     }
}
