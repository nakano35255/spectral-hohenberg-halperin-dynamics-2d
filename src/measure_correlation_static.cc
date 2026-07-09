#include "measure_correlation_static.h"
#include "measure_correlation_static_style.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
namespace {
     constexpr double PI = 3.14159265358979323846;
}
// ---------------------------------------------------------------------- //
CorrelationStaticMeasure::CorrelationStaticMeasure(
     const Params& params,
     const Domain2D& domain,
     std::shared_ptr<const MeasureCommandBase> command
)
     : Measure(params, domain, command),
       spectral_mask_(params, domain),
       local_spectral_size_(domain.spectral_size())
{
     auto static_command = std::dynamic_pointer_cast<const CorrelationStaticMeasureCommand>(command_);
     if (!static_command) {
          throw std::runtime_error("correlation/static measure: invalid command object.");
     }

     nevery_ = static_command->nevery;
     nblock_ = static_command->nblock;
     file_ = static_command->file;
     average_mode_ = static_command->average_mode;
     output_mode_ = static_command->output_mode;
     targets_ = static_command->targets;
     cross_ = static_command->cross;

     nfields_ = static_cast<int>(targets_.size());
     npairs_ = cross_ ? nfields_ * (nfields_ + 1) / 2 : nfields_;

     mode_block_sum_.assign(static_cast<std::size_t>(npairs_) * local_spectral_size_, Complex(0.0, 0.0));
     mode_running_sum_.assign(static_cast<std::size_t>(npairs_) * local_spectral_size_, Complex(0.0, 0.0));

     dk_ = std::min(2.0 * PI / domain_.lx(), 2.0 * PI / domain_.ly());
     if (dk_ <= 0.0) {
          throw std::runtime_error("correlation/static measure: invalid shell spacing.");
     }

     initialize_shells();
}
// ---------------------------------------------------------------------- //
int CorrelationStaticMeasure::shell_index(double k2) const {
     return static_cast<int>(std::floor(std::sqrt(k2) / dk_ + 0.5));
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::initialize_shells() {
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
          throw std::runtime_error("correlation/static found no nonzero active spectral mode.");
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

     shell_block_sum_.assign(static_cast<std::size_t>(npairs_) * static_cast<std::size_t>(nshells_), Complex(0.0, 0.0));
     shell_running_sum_.assign(static_cast<std::size_t>(npairs_) * static_cast<std::size_t>(nshells_), Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
std::size_t CorrelationStaticMeasure::pair_mode_index(int pair_index, std::size_t mode_index) const {
     return static_cast<std::size_t>(pair_index) * local_spectral_size_ + mode_index;
}
// ---------------------------------------------------------------------- //
std::size_t CorrelationStaticMeasure::pair_shell_index(int pair_index, int shell) const {
     return static_cast<std::size_t>(pair_index) * static_cast<std::size_t>(nshells_)
          + static_cast<std::size_t>(shell);
}
// ---------------------------------------------------------------------- //
const Complex* CorrelationStaticMeasure::target_data(const State& state, const CorrelationStaticTarget& target) const {
     switch (target.kind) {
     case CorrelationStaticTargetKind::Rho:
          return state.rho_hat_data();

     case CorrelationStaticTargetKind::Jx:
          return state.jx_hat_data();

     case CorrelationStaticTargetKind::Jy:
          return state.jy_hat_data();

     case CorrelationStaticTargetKind::Psi:
          return state.psi_hat_data(target.component);
     }

     throw std::runtime_error("correlation/static measure: unknown target kind.");
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::accumulate_sample(const State& state) {
     std::vector<const Complex*> fields(static_cast<std::size_t>(nfields_), nullptr);

     for (int field = 0; field < nfields_; ++field) {
          fields[static_cast<std::size_t>(field)] = target_data(state, targets_[static_cast<std::size_t>(field)]);
     }

     const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
     const double spectral_prefactor = 1.0 / (grid_size * grid_size);

     for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
          if (mode.k2 == 0.0) {
               continue;
          }

          const std::size_t mode_index = mode.index;
          const int shell = shell_index(mode.k2);
          const double weight = (mode.gx == 0) ? 1.0 : 2.0;
          const double shell_factor = spectral_prefactor * weight;

          if (cross_) {
               int pair_index = 0;

               for (int first = 0; first < nfields_; ++first) {
                    for (int second = first; second < nfields_; ++second) {
                         const Complex value =
                              fields[static_cast<std::size_t>(first)][mode_index]
                              * std::conj(fields[static_cast<std::size_t>(second)][mode_index]);

                         mode_block_sum_[pair_mode_index(pair_index, mode_index)] +=
                              spectral_prefactor * value;

                         shell_block_sum_[pair_shell_index(pair_index, shell)] +=
                              Complex(shell_factor * std::real(value), 0.0);

                         ++pair_index;
                    }
               }
          } else {
               for (int field = 0; field < nfields_; ++field) {
                    const double value =
                         std::norm(fields[static_cast<std::size_t>(field)][mode_index]);

                    mode_block_sum_[pair_mode_index(field, mode_index)] +=
                         Complex(spectral_prefactor * value, 0.0);

                    shell_block_sum_[pair_shell_index(field, shell)] +=
                         Complex(shell_factor * value, 0.0);
               }
          }
     }

     ++samples_in_block_;
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::finish_block(int step, double time) {
     if (samples_in_block_ <= 0) {
          throw std::runtime_error("CorrelationStaticMeasure tried to finish an empty block.");
     }

     for (std::size_t i = 0; i < mode_block_sum_.size(); ++i) {
          mode_running_sum_[i] += mode_block_sum_[i];
     }

     for (std::size_t i = 0; i < shell_block_sum_.size(); ++i) {
          shell_running_sum_[i] += shell_block_sum_[i];
     }

     ++completed_blocks_;
     running_samples_ += samples_in_block_;

     const int output_samples = (average_mode_ == CorrelationStaticAverageMode::Block) ? samples_in_block_ : running_samples_;
     const std::vector<Complex>& mode_source = (average_mode_ == CorrelationStaticAverageMode::Block) ? mode_block_sum_ : mode_running_sum_;
     const std::vector<Complex>& shell_source = (average_mode_ == CorrelationStaticAverageMode::Block) ? shell_block_sum_ : shell_running_sum_;

     if (output_mode_ == CorrelationStaticOutputMode::TwoD) {
          write_2d_output(step, time, output_samples, mode_source);
     } else {
          write_shell_output(step, time, output_samples, shell_source);
     }

     std::fill(mode_block_sum_.begin(), mode_block_sum_.end(), Complex(0.0, 0.0));
     std::fill(shell_block_sum_.begin(), shell_block_sum_.end(), Complex(0.0, 0.0));
     samples_in_block_ = 0;
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) {
     (void)fft;
     (void)workspace;
     (void)flux;

     ++block_step_;

     if (block_step_ % nevery_ == 0) {
          accumulate_sample(state);
     }

     if (block_step_ == nblock_) {
          finish_block(step, time);
          block_step_ = 0;
     }
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::finalize() {
     if (block_output_.is_open()) {
          block_output_.close();
     }
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::write_pair_column_names(std::ostream& out, bool complex_columns) const {
     if (cross_) {
          for (int first = 0; first < nfields_; ++first) {
               for (int second = first; second < nfields_; ++second) {
                    const std::string name =
                         targets_[static_cast<std::size_t>(first)].name
                         + "_"
                         + targets_[static_cast<std::size_t>(second)].name;

                    if (complex_columns) {
                         out << ' ' << name << "_re"
                             << ' ' << name << "_im";
                    } else {
                         out << ' ' << name;
                    }
               }
          }
     } else {
          for (int field = 0; field < nfields_; ++field) {
               const std::string name =
                    targets_[static_cast<std::size_t>(field)].name
                    + "_"
                    + targets_[static_cast<std::size_t>(field)].name;

               out << ' ' << name;
          }
     }
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::write_header(std::ostream& out) const {
     out << "# measure correlation/static\n";
     out << "# nevery " << nevery_
         << " nblock " << nblock_
         << " mode " << (output_mode_ == CorrelationStaticOutputMode::TwoD ? "2d" : "shell")
         << " average " << (average_mode_ == CorrelationStaticAverageMode::Block ? "block" : "running")
         << " cross " << (cross_ ? "on" : "off")
         << " normalization fft_grid"
         << "\n";
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::open_block_output_if_needed() {
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
          throw std::runtime_error("CorrelationStaticMeasure: failed to open file: " + file_);
     }
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::write_shell_output(int step, double time, int output_samples, const std::vector<Complex>& source) {
     if (average_mode_ == CorrelationStaticAverageMode::Block) {
          open_block_output_if_needed();
     }

     std::vector<double> local_values(source.size(), 0.0);
     for (std::size_t i = 0; i < source.size(); ++i) {
          local_values[i] = source[i].real();
     }

     std::vector<double> global_values;
     if (domain_.rank() == 0) {
          global_values.assign(source.size(), 0.0);
     }

     MPI_Reduce(local_values.data(), global_values.data(), static_cast<int>(local_values.size()), MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

     int write_ok = 1;

     if (domain_.rank() == 0) {
          std::ofstream running_output;
          std::ostream* out = nullptr;

          if (average_mode_ == CorrelationStaticAverageMode::Block) {
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

               (*out) << "# k count";
               write_pair_column_names(*out, false);
               (*out) << "\n";

               for (int shell = 0; shell < nshells_; ++shell) {
                    const double count = shell_weight_counts_[static_cast<std::size_t>(shell)];
                    if (count == 0.0) {
                         continue;
                    }

                    (*out) << static_cast<double>(shell) * dk_
                           << ' ' << count;

                    for (int pair_index = 0; pair_index < npairs_; ++pair_index) {
                         const double sum = global_values[pair_shell_index(pair_index, shell)];
                         const double value = sum / (static_cast<double>(output_samples) * count);

                         (*out) << ' ' << value;
                    }

                    (*out) << '\n';
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
          throw std::runtime_error("CorrelationStaticMeasure: failed to write file: " + file_);
     }
}
// ---------------------------------------------------------------------- //
void CorrelationStaticMeasure::write_2d_output(int step, double time, int output_samples, const std::vector<Complex>& source) {
     if (average_mode_ == CorrelationStaticAverageMode::Block) {
          open_block_output_if_needed();
     }

     std::vector<double> local_values(2 * source.size(), 0.0);
     for (std::size_t i = 0; i < source.size(); ++i) {
          local_values[2 * i] = source[i].real();
          local_values[2 * i + 1] = source[i].imag();
     }

     const int local_count = static_cast<int>(local_values.size());

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

     MPI_Gatherv(
          local_values.data(),
          local_count,
          MPI_DOUBLE,
          gathered.data(),
          counts.data(),
          displs.data(),
          MPI_DOUBLE,
          0,
          domain_.comm()
     );

     int write_ok = 1;

     if (domain_.rank() == 0) {
          const int nkx = domain_.nx_global() / 2 + 1;
          const int nky = domain_.ny_global();
          const std::size_t global_spectral_size =
               static_cast<std::size_t>(nkx) * static_cast<std::size_t>(nky);

          std::vector<Complex> global_values(
               static_cast<std::size_t>(npairs_) * global_spectral_size,
               Complex(0.0, 0.0)
          );

          for (int rank = 0; rank < domain_.size(); ++rank) {
               const Box2D box = domain_.spectral_box_for_rank(rank);
               const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
               const std::size_t rank_spectral_size =
                    static_cast<std::size_t>(box.size[0]) * static_cast<std::size_t>(box.size[1]);
               const std::size_t rank_offset =
                    static_cast<std::size_t>(displs[static_cast<std::size_t>(rank)]);

               for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
                    for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                         const std::size_t local_index =
                              static_cast<std::size_t>(gy - box.low[1]) * local_nkx
                              + static_cast<std::size_t>(gx - box.low[0]);

                         const std::size_t global_index =
                              static_cast<std::size_t>(gy) * static_cast<std::size_t>(nkx)
                              + static_cast<std::size_t>(gx);

                         for (int pair_index = 0; pair_index < npairs_; ++pair_index) {
                              const std::size_t packed_index =
                                   rank_offset
                                   + 2 * (
                                        static_cast<std::size_t>(pair_index) * rank_spectral_size
                                        + local_index
                                   );

                              global_values[static_cast<std::size_t>(pair_index) * global_spectral_size + global_index] =
                                   Complex(gathered[packed_index], gathered[packed_index + 1])
                                   / static_cast<double>(output_samples);
                         }
                    }
               }
          }

          std::ofstream running_output;
          std::ostream* out = nullptr;

          if (average_mode_ == CorrelationStaticAverageMode::Block) {
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

               (*out) << "# kx ky";
               write_pair_column_names(*out, cross_);
               (*out) << "\n";

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

                         const std::size_t global_index =
                              static_cast<std::size_t>(gy) * static_cast<std::size_t>(nkx)
                              + static_cast<std::size_t>(gx);

                         (*out) << kx << ' ' << ky;

                         for (int pair_index = 0; pair_index < npairs_; ++pair_index) {
                              const Complex value =
                                   global_values[static_cast<std::size_t>(pair_index) * global_spectral_size + global_index];

                              if (cross_) {
                                   (*out) << ' ' << value.real()
                                          << ' ' << value.imag();
                              } else {
                                   (*out) << ' ' << value.real();
                              }
                         }

                         (*out) << '\n';
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
          throw std::runtime_error("CorrelationStaticMeasure: failed to write file: " + file_);
     }
}
// ---------------------------------------------------------------------- //

