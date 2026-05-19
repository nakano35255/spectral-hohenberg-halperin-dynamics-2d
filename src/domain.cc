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
    bool is_power_of_two(int value) {
        return value > 0 && (value & (value - 1)) == 0;
    }
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
void Domain2D::validate_grid(int gnx, int gny, double glx, double gly, int size) {
    if (gnx <= 0 || gny <= 0) {
        throw std::runtime_error("Domain2D requires positive grid sizes.");
    }
    if (!is_power_of_two(gnx) || !is_power_of_two(gny)) {
        throw std::runtime_error("Domain2D requires power-of-two grid sizes.");
    }
    if (glx <= 0.0 || gly <= 0.0) {
        throw std::runtime_error("Domain2D requires positive domain lengths.");
    }
    if (size <= 0) {
        throw std::runtime_error("Domain2D received an invalid MPI size.");
    }
    if (!is_power_of_two(size)) {
        throw std::runtime_error("Domain2D requires a power-of-two MPI size.");
    }
    if (size > gny) {
        throw std::runtime_error("Domain2D y-slab decomposition requires size <= Ny.");
    }
    if (gny % size != 0) {
        throw std::runtime_error("Domain2D y-slab decomposition requires Ny to be divisible by MPI size.");
    }
    if (gny / size < 2) {
        throw std::runtime_error("Domain2D y-slab decomposition requires at least 2 local y grid points per rank.");
    }
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::make_physical_y_slab_box(
    int gnx, int gny, double glx, double gly, int rank, int size
) {
    validate_grid(gnx, gny, glx, gly, size);

    const int local_ny = gny / size;
    const int y0 = rank * local_ny;
    const int y1 = y0 + local_ny - 1;

    return {{0, y0, 0}, {gnx - 1, y1, 0}};
}
// ---------------------------------------------------------------------- //
Box2D Domain2D::make_r2c_spectral_box(const Box2D& physical_box) {
    return physical_box.r2c(R2C_DIRECTION);
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
      physical_box_(make_physical_y_slab_box(gnx, gny, glx, gly, rank_, size_)),
      spectral_box_(make_r2c_spectral_box(physical_box_)) {}
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
