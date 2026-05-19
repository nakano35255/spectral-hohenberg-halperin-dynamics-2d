#include "measure_snapshot.h"
#include "measure_snapshot_style.h"

#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
SnapshotMeasure::SnapshotMeasure(const Params& params, std::shared_ptr<const MeasureCommandBase> command) :
       Measure(params, command),
       num_components_(params.physics.num_components),
       step_width_(compute_step_width(params.total_run_steps()))
{
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
     out << "# nx " << nx << " ny " << ny << " components " << num_components_ << "\n";
     out << "# columns x y field component value\n";
     out << std::scientific << std::setprecision(16);

     for (int gy = 0; gy < ny; ++gy) {
          for (int gx = 0; gx < nx; ++gx) {
               const int grid_index = gy * nx + gx;
               for (int component = 0; component < num_components_; ++component) {
                    out << gx << ' ' << gy << " rho " << component << ' '
                        << global_data[component * field_size + grid_index] << '\n';
               }
               out << gx << ' ' << gy << " jx -1 "
                   << global_data[num_components_ * field_size + grid_index] << '\n';
               out << gx << ' ' << gy << " jy -1 "
                   << global_data[(num_components_ + 1) * field_size + grid_index] << '\n';
          }
     }
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::write_spectral_snapshot(
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

     const int nkx = domain.nx_global() / 2 + 1;
     const int nky = domain.ny_global();
     const int field_size = nkx * nky;

     out << "# snapshot space spectral\n";
     out << "# step " << step << " time " << std::scientific << std::setprecision(16) << time << "\n";
     out << "# nkx " << nkx << " nky " << nky << " components " << num_components_ << "\n";
     out << "# columns kx_index ky_index field component real imag\n";
     out << std::scientific << std::setprecision(16);

     for (int ky = 0; ky < nky; ++ky) {
          for (int kx = 0; kx < nkx; ++kx) {
               const int grid_index = ky * nkx + kx;
               for (int component = 0; component < num_components_; ++component) {
                    const int offset = 2 * (component * field_size + grid_index);
                    out << kx << ' ' << ky << " rho " << component << ' '
                        << global_data[offset] << ' ' << global_data[offset + 1] << '\n';
               }

               const int jx_offset = 2 * (num_components_ * field_size + grid_index);
               const int jy_offset = 2 * ((num_components_ + 1) * field_size + grid_index);
               out << kx << ' ' << ky << " jx -1 "
                   << global_data[jx_offset] << ' ' << global_data[jx_offset + 1] << '\n';
               out << kx << ' ' << ky << " jy -1 "
                   << global_data[jy_offset] << ' ' << global_data[jy_offset + 1] << '\n';
          }
     }
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::observe_physical(
     const State& state,
     PhysicalStateBuffer& physical,
     FourierTransform2D& fft,
     const Domain2D& domain,
     int step,
     double time
) const {
     const int nfields = num_components_ + 2;
     const int nx = domain.nx_global();
     const int local_ny = domain.physical_box().size[1];
     const int local_count = nfields * nx * local_ny;

     fft.backward_many(nfields, state.data(), physical.data());

     std::vector<double> gathered;
     if (domain.rank() == 0) {
          gathered.resize(static_cast<std::size_t>(local_count * domain.size()));
     }

     MPI_Gather(
          physical.data(),
          local_count,
          MPI_DOUBLE,
          gathered.data(),
          local_count,
          MPI_DOUBLE,
          0,
          domain.comm()
     );

     if (domain.rank() != 0) {
          return;
     }

     std::vector<double> global_data(static_cast<std::size_t>(nfields * nx * domain.ny_global()));
     const int local_field_size = nx * local_ny;
     const int global_field_size = nx * domain.ny_global();

     for (int rank = 0; rank < domain.size(); ++rank) {
          const int y0 = rank * local_ny;
          for (int field = 0; field < nfields; ++field) {
               for (int ly = 0; ly < local_ny; ++ly) {
                    for (int gx = 0; gx < nx; ++gx) {
                         const int local_index =
                              rank * local_count + field * local_field_size + ly * nx + gx;
                         const int global_index =
                              field * global_field_size + (y0 + ly) * nx + gx;
                         global_data[global_index] = gathered[local_index];
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
     const int nfields = num_components_ + 2;
     const int nkx = domain.nx_global() / 2 + 1;
     const int local_ny = domain.spectral_box().size[1];
     const int local_complex_count = nfields * nkx * local_ny;

     std::vector<double> local_data(static_cast<std::size_t>(2 * local_complex_count));
     const Complex* spectral = state.data();
     for (int i = 0; i < local_complex_count; ++i) {
          local_data[2 * i] = spectral[i].real();
          local_data[2 * i + 1] = spectral[i].imag();
     }

     std::vector<double> gathered;
     if (domain.rank() == 0) {
          gathered.resize(static_cast<std::size_t>(2 * local_complex_count * domain.size()));
     }

     MPI_Gather(
          local_data.data(),
          2 * local_complex_count,
          MPI_DOUBLE,
          gathered.data(),
          2 * local_complex_count,
          MPI_DOUBLE,
          0,
          domain.comm()
     );

     if (domain.rank() != 0) {
          return;
     }

     std::vector<double> global_data(static_cast<std::size_t>(2 * nfields * nkx * domain.ny_global()));
     const int local_field_size = nkx * local_ny;
     const int global_field_size = nkx * domain.ny_global();
     const int local_double_count = 2 * local_complex_count;

     for (int rank = 0; rank < domain.size(); ++rank) {
          const int y0 = rank * local_ny;
          for (int field = 0; field < nfields; ++field) {
               for (int ly = 0; ly < local_ny; ++ly) {
                    for (int kx = 0; kx < nkx; ++kx) {
                         const int local_complex_index =
                              field * local_field_size + ly * nkx + kx;
                         const int global_complex_index =
                              field * global_field_size + (y0 + ly) * nkx + kx;
                         global_data[2 * global_complex_index] =
                              gathered[rank * local_double_count + 2 * local_complex_index];
                         global_data[2 * global_complex_index + 1] =
                              gathered[rank * local_double_count + 2 * local_complex_index + 1];
                    }
               }
          }
     }

     write_spectral_snapshot(snapshot_filename(step), global_data, domain, step, time);
}
// ---------------------------------------------------------------------- //
void SnapshotMeasure::observe(const State& state, PhysicalStateBuffer& physical, FourierTransform2D& fft, const Domain2D& domain, int step, double time) {
     if (nevery_ <= 0 || step % nevery_ != 0) {
          return;
     }

     if (space_ == "physical") {
          observe_physical(state, physical, fft, domain, step, time);
     }
     else if (space_ == "spectral") {
          observe_spectral(state, domain, step, time);
     }
     else {
          throw std::runtime_error("Unknown snapshot space: " + space_);
     }
}
// ---------------------------------------------------------------------- //
