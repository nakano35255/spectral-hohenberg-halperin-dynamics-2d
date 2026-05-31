#include "measure_snapshot.h"
#include "measure_snapshot_style.h"

#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
SnapshotMeasure::SnapshotMeasure(
    const Params& params,
    const Domain2D& domain,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    num_order_parameters_(params.physics.num_order_parameters),
    num_fields_(params.physics.num_order_parameters + 3),
    step_width_(compute_step_width(params.total_run_steps()))
{
     if (num_order_parameters_ < 0) {
          throw std::runtime_error("SnapshotMeasure requires a nonnegative number of order parameters.");
     }

     auto cfg = std::dynamic_pointer_cast<const SnapshotMeasureCommand>(command);
     if (!cfg) {
          throw std::runtime_error("SnapshotMeasure: invalid command type.");
     }

     nevery_ = cfg->nevery;
     file_ = cfg->file;
     space_ = cfg->space;
}
// ---------------------------------------------------------------------- //
int SnapshotMeasure::compute_step_width(int max_steps) {
     int digits = 1;
     while (max_steps >= 10) {
          max_steps /= 10;
          ++digits;
     }
     return digits + 1;
}
// ---------------------------------------------------------------------- //
std::string SnapshotMeasure::measurement_step_string(int value, int width) {
     std::ostringstream oss;
     oss << std::setw(width) << std::setfill('0') << value;
     return oss.str();
}
// ---------------------------------------------------------------------- //
std::string SnapshotMeasure::snapshot_filename(int step) const {
     return file_ + "_step" + measurement_step_string(step, step_width_) + ".snapshot";
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::write_physical_snapshot(
     const std::string& filename,
     const std::vector<double>& global_data,
     const Domain2D& domain,
     int step,
     double time
) const {
     std::ofstream out(filename);
     if (!out) {
          throw std::runtime_error("SnapshotMeasure: cannot open file: " + filename);
     }

     const int nx = domain.nx_global();
     const int ny = domain.ny_global();
     const int field_size = nx * ny;

     out << "# snapshot space physical\n";
     out << "# step " << step << " time " << std::scientific << std::setprecision(16) << time << "\n";
     out << "# nx " << nx << " ny " << ny << " order_parameters " << num_order_parameters_ << "\n";
     out << "# x y rho";
     for (int component = 0; component < num_order_parameters_; ++component) {
          out << " psi_com" << component;
     }
     out << " jx jy\n";
     out << std::scientific << std::setprecision(16);

     for (int gy = 0; gy < ny; ++gy) {
          for (int gx = 0; gx < nx; ++gx) {
               const int grid_index = gy * nx + gx;
               out << gx << ' ' << gy;
               out << ' ' << global_data[grid_index];
               for (int component = 0; component < num_order_parameters_; ++component) {
                    out << ' ' << global_data[(component + 1) * field_size + grid_index];
               }
               out << ' ' << global_data[(num_order_parameters_ + 1) * field_size + grid_index]
                   << ' ' << global_data[(num_order_parameters_ + 2) * field_size + grid_index]
                   << '\n';
          }
     }
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::write_spectral_snapshot(const std::string& filename, const std::vector<double>& global_data, const Domain2D& domain, int step, double time) const {
     std::ofstream out(filename);
     if (!out) {
          throw std::runtime_error("SnapshotMeasure: cannot open file: " + filename);
     }

     const int nkx = domain.nx_global() / 2 + 1;
     const int nky = domain.ny_global();
     const int field_size = nkx * nky;

     out << "# snapshot space spectral\n";
     out << "# step " << step << " time " << std::scientific << std::setprecision(16) << time << "\n";
     out << "# nkx " << nkx << " nky " << nky << " order_parameters " << num_order_parameters_ << "\n";
     out << "# kx_index ky_index rho_real rho_imag";
     for (int component = 0; component < num_order_parameters_; ++component) {
          out << " psi_com" << component << "_real psi_com" << component << "_imag";
     }
     out << " jx_real jx_imag jy_real jy_imag\n";
     out << std::scientific << std::setprecision(16);

     for (int ky = 0; ky < nky; ++ky) {
          for (int kx = 0; kx < nkx; ++kx) {
               const int grid_index = ky * nkx + kx;
               out << kx << ' ' << ky;
               const int rho_offset = 2 * grid_index;
               out << ' ' << global_data[rho_offset]
                   << ' ' << global_data[rho_offset + 1];

               for (int component = 0; component < num_order_parameters_; ++component) {
                    const int offset = 2 * ((component + 1) * field_size + grid_index);
                    out << ' ' << global_data[offset]
                        << ' ' << global_data[offset + 1];
               }

               const int jx_offset = 2 * ((num_order_parameters_ + 1) * field_size + grid_index);
               const int jy_offset = 2 * ((num_order_parameters_ + 2) * field_size + grid_index);
               out << ' ' << global_data[jx_offset]
                   << ' ' << global_data[jx_offset + 1]
                   << ' ' << global_data[jy_offset]
                   << ' ' << global_data[jy_offset + 1]
                   << '\n';
          }
     }
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::observe_physical(
     PhysicalStateBuffer& physical,
     const Domain2D& domain,
     int step,
     double time
) const {
     const int nfields = num_fields_;
     const int nx = domain.nx_global();
     const int ny = domain.ny_global();

     const Box2D& local_box = domain.physical_box();
     const int local_count = nfields * local_box.size[0] * local_box.size[1];

     std::vector<int> counts;
     if (domain.rank() == 0) {
          counts.resize(static_cast<std::size_t>(domain.size()));
     }

     MPI_Gather(
          &local_count, 1, MPI_INT,
          counts.data(), 1, MPI_INT,
          0, domain.comm()
     );

     std::vector<int> displs;
     std::vector<double> gathered;

     if (domain.rank() == 0) {
          displs.resize(static_cast<std::size_t>(domain.size()), 0);
          int total = 0;
          for (int rank = 0; rank < domain.size(); ++rank) {
               displs[rank] = total;
               total += counts[rank];
          }
          gathered.resize(static_cast<std::size_t>(total));
     }

     MPI_Gatherv(
          physical.data(), local_count, MPI_DOUBLE,
          gathered.data(), counts.data(), displs.data(), MPI_DOUBLE,
          0, domain.comm()
     );

     if (domain.rank() != 0) {
          return;
     }

     std::vector<double> global_data(static_cast<std::size_t>(nfields * nx * ny));
     const int global_field_size = nx * ny;

     for (int rank = 0; rank < domain.size(); ++rank) {
          const Box2D box = domain.physical_box_for_rank(rank);
          const int local_nx = box.size[0];
          const int local_ny = box.size[1];
          const int local_field_size = local_nx * local_ny;
          const int rank_offset = displs[rank];

          for (int field = 0; field < nfields; ++field) {
               for (int ly = 0; ly < local_ny; ++ly) {
                    const int gy = box.low[1] + ly;

                    for (int lx = 0; lx < local_nx; ++lx) {
                         const int gx = box.low[0] + lx;

                         const int local_index =
                              field * local_field_size + ly * local_nx + lx;
                         const int global_index =
                              field * global_field_size + gy * nx + gx;

                         global_data[global_index] =
                              gathered[rank_offset + local_index];
                    }
               }
          }
     }

     write_physical_snapshot(snapshot_filename(step), global_data, domain, step, time);
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::observe_spectral(
     const State& state,
     const Domain2D& domain,
     int step,
     double time
) const {
     const int nfields = num_fields_;
     const int nkx = domain.nx_global() / 2 + 1;
     const int nky = domain.ny_global();

     const Box2D& local_box = domain.spectral_box();
     const int local_complex_count =
          nfields * local_box.size[0] * local_box.size[1];

     std::vector<double> local_data(static_cast<std::size_t>(2 * local_complex_count));

     const Complex* spectral = state.data();
     for (int i = 0; i < local_complex_count; ++i) {
          local_data[2 * i] = spectral[i].real();
          local_data[2 * i + 1] = spectral[i].imag();
     }

     const int local_double_count = 2 * local_complex_count;

     std::vector<int> counts;
     if (domain.rank() == 0) {
          counts.resize(static_cast<std::size_t>(domain.size()));
     }

     MPI_Gather(
          &local_double_count, 1, MPI_INT,
          counts.data(), 1, MPI_INT,
          0, domain.comm()
     );

     std::vector<int> displs;
     std::vector<double> gathered;

     if (domain.rank() == 0) {
          displs.resize(static_cast<std::size_t>(domain.size()), 0);
          int total = 0;
          for (int rank = 0; rank < domain.size(); ++rank) {
               displs[rank] = total;
               total += counts[rank];
          }
          gathered.resize(static_cast<std::size_t>(total));
     }

     MPI_Gatherv(
          local_data.data(), local_double_count, MPI_DOUBLE,
          gathered.data(), counts.data(), displs.data(), MPI_DOUBLE,
          0, domain.comm()
     );

     if (domain.rank() != 0) {
          return;
     }

     std::vector<double> global_data(
          static_cast<std::size_t>(2 * nfields * nkx * nky)
     );

     const int global_field_size = nkx * nky;

     for (int rank = 0; rank < domain.size(); ++rank) {
          const Box2D box = domain.spectral_box_for_rank(rank);
          const int local_nkx = box.size[0];
          const int local_nky = box.size[1];
          const int local_field_size = local_nkx * local_nky;
          const int rank_offset = displs[rank];

          for (int field = 0; field < nfields; ++field) {
               for (int ly = 0; ly < local_nky; ++ly) {
                    const int ky = box.low[1] + ly;

                    for (int lx = 0; lx < local_nkx; ++lx) {
                         const int kx = box.low[0] + lx;

                         const int local_complex_index =
                              field * local_field_size + ly * local_nkx + lx;
                         const int global_complex_index =
                              field * global_field_size + ky * nkx + kx;

                         global_data[2 * global_complex_index] =
                              gathered[rank_offset + 2 * local_complex_index];
                         global_data[2 * global_complex_index + 1] =
                              gathered[rank_offset + 2 * local_complex_index + 1];
                    }
               }
          }
     }

     write_spectral_snapshot(snapshot_filename(step), global_data, domain, step, time);
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::observe(
     const State& state,
     FourierTransform2D& fft,
     MeasureWorkspace& workspace,
     const FluxBuffer& /*flux*/,
     int step,
     double time
) {
     if (nevery_ <= 0 || step % nevery_ != 0) {
          return;
     }

     if (space_ == "physical") {
          PhysicalStateBuffer& physical = workspace.ensure_physical(state, fft);
          observe_physical(physical, domain_, step, time);
     }
     else if (space_ == "spectral") {
          observe_spectral(state, domain_, step, time);
     }
     else {
          throw std::runtime_error("Unknown snapshot space: " + space_);
     }
}
// ---------------------------------------------------------------------- //
