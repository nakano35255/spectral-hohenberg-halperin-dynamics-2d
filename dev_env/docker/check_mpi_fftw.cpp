#define OMPI_SKIP_MPICXX 1

#include <mpi.h>
#include <fftw3-mpi.h>

#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const ptrdiff_t n0 = 8;
    const ptrdiff_t n1 = 8;
    ptrdiff_t local_n0 = 0;
    ptrdiff_t local_0_start = 0;

    const ptrdiff_t alloc_local = fftw_mpi_local_size_2d(
        n0, n1, MPI_COMM_WORLD, &local_n0, &local_0_start);

    fftw_complex* data = fftw_alloc_complex(alloc_local);
    if (data == nullptr) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        throw std::runtime_error("fftw_alloc_complex failed");
    }

    for (ptrdiff_t i = 0; i < alloc_local; ++i) {
        data[i][0] = 0.0;
        data[i][1] = 0.0;
    }

    fftw_plan plan = fftw_mpi_plan_dft_2d(
        n0, n1, data, data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);

    if (plan == nullptr) {
        fftw_free(data);
        MPI_Abort(MPI_COMM_WORLD, 2);
        throw std::runtime_error("fftw_mpi_plan_dft_2d failed");
    }

    fftw_execute(plan);

    if (rank == 0) {
        std::cout << "FFTW MPI check OK with " << size << " ranks\n";
    }

    fftw_destroy_plan(plan);
    fftw_free(data);
    MPI_Finalize();
    return 0;
}
