#ifndef SFI_DOMAIN_H
#define SFI_DOMAIN_H

#define OMPI_SKIP_MPICXX 1
#include "simulationinfo.h"

#include <mpi.h>
#include <heffte.h>
#include <cstddef>


using Box2D = heffte::box3d<>;

class Domain2D {
public:
    static constexpr int R2C_DIRECTION = 0;

private:
    MPI_Comm comm_ = MPI_COMM_WORLD;
    int rank_ = 0;
    int size_ = 1;

    int global_nx_ = 0;
    int global_ny_ = 0;
    double global_lx_ = 0.0;
    double global_ly_ = 0.0;

    Box2D physical_box_;
    Box2D spectral_box_;

    static int mpi_rank(MPI_Comm comm);
    static int mpi_size(MPI_Comm comm);
    static void validate_grid(int gnx, int gny, double glx, double gly, int size);
    static Box2D make_physical_y_slab_box(
        int gnx, int gny, double glx, double gly, int rank, int size
    );
    static Box2D make_r2c_spectral_box(const Box2D& physical_box);

public:
    Domain2D(int gnx, int gny, double glx, double gly, MPI_Comm comm = MPI_COMM_WORLD);
    Domain2D(const Params& params, MPI_Comm comm = MPI_COMM_WORLD);

    MPI_Comm comm() const { return comm_; }
    int rank() const { return rank_; }
    int size() const { return size_; }

    int nx_global() const { return global_nx_; }
    int ny_global() const { return global_ny_; }

    double lx() const { return global_lx_; }
    double ly() const { return global_ly_; }

    const Box2D& physical_box() const { return physical_box_; }
    const Box2D& spectral_box() const { return spectral_box_; }

    std::size_t physical_size() const;
    std::size_t spectral_size() const;
    std::size_t physical_local_index(int gx, int gy) const;
    std::size_t spectral_local_index(int gx, int gy) const;

    double kx(int gx) const;
    double ky(int gy) const;
};

#endif
