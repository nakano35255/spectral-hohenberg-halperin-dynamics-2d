#define OMPI_SKIP_MPICXX 1

#include <mpi.h>
#include <heffte.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

namespace {

using complex_t = std::complex<double>;

int global_index(int i0, int i1, int i2) {
    const int n0 = 4;
    const int n1 = 4;
    return i2 * n0 * n1 + i1 * n0 + i0;
}

int local_index(const heffte::box3d<>& box, int i0, int i1, int i2) {
    const int local_plane = box.size[0] * box.size[1];
    const int local_stride = box.size[0];
    return (i2 - box.low[2]) * local_plane
         + (i1 - box.low[1]) * local_stride
         + (i0 - box.low[0]);
}

complex_t initial_value(int i0, int i1, int i2) {
    const double id = static_cast<double>(global_index(i0, i1, i2));
    return complex_t(id, -0.25 * id);
}

double local_max_error(const heffte::box3d<>& box,
                       const std::vector<complex_t>& values) {
    double err = 0.0;

    for (int i2 = box.low[2]; i2 <= box.high[2]; ++i2) {
        for (int i1 = box.low[1]; i1 <= box.high[1]; ++i1) {
            for (int i0 = box.low[0]; i0 <= box.high[0]; ++i0) {
                const int loc = local_index(box, i0, i1, i2);
                err = std::max(err, std::abs(values[loc] - initial_value(i0, i1, i2)));
            }
        }
    }

    return err;
}

} // namespace

int run_check(MPI_Comm comm) {
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (size != 2) {
        if (rank == 0) {
            std::cerr << "This heFFTe check expects exactly 2 MPI ranks.\n";
            std::cerr << "Run: mpirun --allow-run-as-root -np 2 ./out.exe\n";
        }
        return 1;
    }

    const heffte::box3d<> left_box = {{0, 0, 0}, {3, 3, 1}};
    const heffte::box3d<> right_box = {{0, 0, 2}, {3, 3, 3}};
    const heffte::box3d<> my_box = (rank == 0) ? left_box : right_box;

    heffte::fft3d<heffte::backend::fftw> fft(my_box, my_box, comm);

    std::vector<complex_t> input(fft.size_inbox());
    std::vector<complex_t> spectral(fft.size_outbox());
    std::vector<complex_t> recovered(fft.size_inbox());

    for (int i2 = my_box.low[2]; i2 <= my_box.high[2]; ++i2) {
        for (int i1 = my_box.low[1]; i1 <= my_box.high[1]; ++i1) {
            for (int i0 = my_box.low[0]; i0 <= my_box.high[0]; ++i0) {
                input[local_index(my_box, i0, i1, i2)] = initial_value(i0, i1, i2);
            }
        }
    }

    fft.forward(input.data(), spectral.data());
    fft.backward(spectral.data(), recovered.data());

    const double grid_size = 4.0 * 4.0 * 4.0;
    for (auto& value : recovered) {
        value /= grid_size;
    }

    const double local_error = local_max_error(my_box, recovered);
    double global_error = 0.0;
    MPI_Reduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    int status = 0;
    if (rank == 0) {
        std::cout << std::scientific;
        std::cout << "heFFTe FFTW max error: " << global_error << '\n';

        if (global_error < 1.0e-10) {
            std::cout << "heFFTe FFTW check OK with " << size << " ranks\n";
        } else {
            std::cerr << "heFFTe FFTW check failed\n";
            status = 1;
        }
    }

    MPI_Bcast(&status, 1, MPI_INT, 0, comm);
    return status;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    const int status = run_check(MPI_COMM_WORLD);

    MPI_Finalize();
    return status;
}
