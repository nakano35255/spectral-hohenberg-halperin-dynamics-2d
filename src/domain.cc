#include "domain.h"

#include "simulationinfo.h"

#include <cmath>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
namespace {
    // ---------------------------------------------------------------------- //
    constexpr double PI = 3.14159265358979323846;
    // ---------------------------------------------------------------------- //
    std::size_t box_size_2d(const Box2D& box) {
        if (box.empty()) {
            return 0;
        }
        return static_cast<std::size_t>(box.size[0]) * static_cast<std::size_t>(box.size[1]);
    }
    // ---------------------------------------------------------------------- //
    std::size_t box_local_index_2d(const Box2D& box, int gx, int gy) {
        if (box.empty() || gx < box.low[0] || gx > box.high[0] || gy < box.low[1] || gy > box.high[1]) {
            throw std::out_of_range("Box2D index out of range: (" + std::to_string(gx) + ", " + std::to_string(gy) + ")");
        }

        const int local_x = gx - box.low[0];
        const int local_y = gy - box.low[1];

        return static_cast<std::size_t>(local_y) * static_cast<std::size_t>(box.size[0]) + static_cast<std::size_t>(local_x);
    }
    // ---------------------------------------------------------------------- //
}
// ---------------------------------------------------------------------- //
int Domain2D::mpi_rank(MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank;
}
// ---------------------------------------------------------------------- //
int Domain2D::mpi_size(MPI_Comm comm) {
    int size = 1;
    MPI_Comm_size(comm, &size);
    return size;
}
// ---------------------------------------------------------------------- //
void Domain2D::validate_grid(int gnx, int gny, int size) {
    if (gnx <= 0 || gny <= 0) {
        throw std::runtime_error("Domain2D requires positive grid sizes.");
    }
    if (size <= 0) {
        throw std::runtime_error("Domain2D received an invalid MPI size.");
    }

    const int spectral_nx = gnx / 2 + 1;

    if (size > gnx) {
        throw std::runtime_error("Domain2D x-slab decomposition requires size <= Nx.");
    }
    if (size > spectral_nx) {
        throw std::runtime_error("Domain2D spectral x-slab decomposition requires size <= Nx/2+1.");
    }
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::make_x_slab_box(int nx, int ny, int rank, int size) {
    if (nx <= 0 || ny <= 0) {
        throw std::runtime_error("x-slab box requires positive sizes.");
    }
    if (rank < 0 || rank >= size) {
        throw std::runtime_error("x-slab box received invalid rank.");
    }
    if (size > nx) {
        throw std::runtime_error("x-slab box requires size <= nx.");
    }

    const int base = nx / size;
    const int remainder = nx % size;

    const int local_nx = base + (rank < remainder ? 1 : 0);
    const int x0 = rank * base + std::min(rank, remainder);
    const int x1 = x0 + local_nx - 1;

    return {{x0, 0, 0}, {x1, ny - 1, 0}};
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::make_physical_x_slab_box(int gnx, int gny, int rank, int size) {
    validate_grid(gnx, gny, size);
    return make_x_slab_box(gnx, gny, rank, size);
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::make_r2c_spectral_x_slab_box(int gnx, int gny, int rank, int size) {
    const int spectral_nx = gnx / 2 + 1;
    return make_x_slab_box(spectral_nx, gny, rank, size);
}
// ---------------------------------------------------------------------- //
Domain2D::Domain2D(int gnx, int gny, double glx, double gly, MPI_Comm comm)
    : comm_(comm),
      rank_(mpi_rank(comm)),
      size_(mpi_size(comm)),
      global_nx_(gnx),
      global_ny_(gny),
      global_lx_(glx),
      global_ly_(gly),
      physical_box_(make_physical_x_slab_box(gnx, gny, rank_, size_)),
      spectral_box_(make_r2c_spectral_x_slab_box(gnx, gny, rank_, size_)) {}
// ---------------------------------------------------------------------- //
Domain2D::Domain2D(const Params& params, MPI_Comm comm)
    : Domain2D(
          params.grid.num_grid[0],
          params.grid.num_grid[1],
          params.grid.length[0],
          params.grid.length[1],
          comm
      ) {}
// ---------------------------------------------------------------------- //
std::size_t Domain2D::physical_size() const {
    return box_size_2d(physical_box_);
}
// ---------------------------------------------------------------------- //
std::size_t Domain2D::spectral_size() const {
    return box_size_2d(spectral_box_);
}
// ---------------------------------------------------------------------- //
std::size_t Domain2D::physical_local_index(int gx, int gy) const {
    return box_local_index_2d(physical_box_, gx, gy);
}
// ---------------------------------------------------------------------- //
std::size_t Domain2D::spectral_local_index(int gx, int gy) const {
    return box_local_index_2d(spectral_box_, gx, gy);
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::physical_box_for_rank(int rank) const {
    return make_physical_x_slab_box(global_nx_, global_ny_, rank, size_);
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::spectral_box_for_rank(int rank) const {
    return make_r2c_spectral_x_slab_box(global_nx_, global_ny_, rank, size_);
}
// ---------------------------------------------------------------------- //
double Domain2D::kx(int gx) const {
    if (gx < 0 || gx > global_nx_ / 2) {
        throw std::out_of_range("R2C kx index out of range: " + std::to_string(gx));
    }
    return 2.0 * PI * static_cast<double>(gx) / global_lx_;
}
// ---------------------------------------------------------------------- //
double Domain2D::ky(int gy) const {
    const int wrapped = (gy <= global_ny_ / 2) ? gy : gy - global_ny_;
    return 2.0 * PI * static_cast<double>(wrapped) / global_ly_;
}
// ---------------------------------------------------------------------- //
