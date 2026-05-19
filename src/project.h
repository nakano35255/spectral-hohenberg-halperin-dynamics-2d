#ifndef PROJECT_H
#define PROJECT_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "buffer_physical_state.h"
#include "fourier_transform.h"
#include "monitor.h"
#include "measure_registry_builtin.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class Project {
private:
    Params params;

    Domain2D domain;
    State state;
    PhysicalStateBuffer buf_physical_state;
    FourierTransform2D fourier;

    MeasureRegistry registry;
    SimulationMonitor monitor;
    int step;
    int run_index;
    double time;

    void run_fft_sine_mode_test() {
        const int nx = domain.nx_global();
        const int ny = domain.ny_global();
        const double two_pi = 6.28318530717958647692;

        const Box2D& physical_box = domain.physical_box();
        for (int gy = physical_box.low[1]; gy <= physical_box.high[1]; ++gy) {
            const double y = domain.ly() * static_cast<double>(gy) / static_cast<double>(ny);
            for (int gx = physical_box.low[0]; gx <= physical_box.high[0]; ++gx) {
                const double x = domain.lx() * static_cast<double>(gx) / static_cast<double>(nx);
                buf_physical_state.rho(0, gx, gy) = std::sin(two_pi * x / domain.lx())
                                                   * std::sin(2.0 * two_pi * y / domain.ly());
            }
        }

        fourier.forward(buf_physical_state.rho_data(0), state.rho_hat_data(0));

        const double expected_amplitude = 0.25 * static_cast<double>(nx) * static_cast<double>(ny);
        double local_max_error = 0.0;

        const Box2D& spectral_box = domain.spectral_box();
        for (int gy = spectral_box.low[1]; gy <= spectral_box.high[1]; ++gy) {
            for (int gx = spectral_box.low[0]; gx <= spectral_box.high[0]; ++gx) {
                Complex expected(0.0, 0.0);
                if (gx == 1 && gy == 2) {
                    expected = Complex(-expected_amplitude, 0.0);
                }
                else if (gx == 1 && gy == ny - 2) {
                    expected = Complex(expected_amplitude, 0.0);
                }

                const Complex actual = state.rho_hat(0, gx, gy);
                local_max_error = std::max(local_max_error, std::abs(actual - expected));
            }
        }

        double global_max_error = 0.0;
        MPI_Allreduce(
            &local_max_error,
            &global_max_error,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            domain.comm()
        );

        if (global_max_error > 1.0e-10) {
            throw std::runtime_error("FFT sine mode test failed.");
        }

        if (nx <= 8 && ny <= 8) {
            std::vector<double> local_modes(2 * domain.spectral_size(), 0.0);
            std::size_t local_mode_index = 0;
            for (int gy = spectral_box.low[1]; gy <= spectral_box.high[1]; ++gy) {
                for (int gx = spectral_box.low[0]; gx <= spectral_box.high[0]; ++gx) {
                    const Complex value = state.rho_hat(0, gx, gy);
                    local_modes[2 * local_mode_index] = value.real();
                    local_modes[2 * local_mode_index + 1] = value.imag();
                    ++local_mode_index;
                }
            }

            const int local_count = static_cast<int>(local_modes.size());
            std::vector<double> all_modes;
            if (domain.rank() == 0) {
                all_modes.resize(static_cast<std::size_t>(local_count) * static_cast<std::size_t>(domain.size()));
            }

            MPI_Gather(
                local_modes.data(),
                local_count,
                MPI_DOUBLE,
                domain.rank() == 0 ? all_modes.data() : nullptr,
                local_count,
                MPI_DOUBLE,
                0,
                domain.comm()
            );

            if (domain.rank() == 0) {
                std::ostringstream os;
                os << "\n[FFT sine mode test]\n";
                os << "  input: sin(2*pi*x/Lx) * sin(2*2*pi*y/Ly)\n";
                os << "  expected nonzero R2C modes: (kx,ky)=(1,2), (1," << ny - 2 << ")\n";
                os << "  max spectral error: " << std::scientific << global_max_error << "\n\n";

                const int local_ny = ny / domain.size();
                const int spectral_nx = nx / 2 + 1;
                for (int rank = 0; rank < domain.size(); ++rank) {
                    const int y0 = rank * local_ny;
                    const int y1 = y0 + local_ny - 1;
                    std::size_t offset = static_cast<std::size_t>(rank)
                                       * static_cast<std::size_t>(local_count);

                    os << "Rank " << rank << " spectral box: "
                       << "x=[0," << nx / 2 << "] "
                       << "y=[" << y0 << "," << y1 << "]\n";

                    std::size_t mode_index = 0;
                    for (int gy = y0; gy <= y1; ++gy) {
                        for (int gx = 0; gx < spectral_nx; ++gx) {
                            const double real = all_modes[offset + 2 * mode_index];
                            const double imag = all_modes[offset + 2 * mode_index + 1];
                            os << "  (" << std::setw(1) << gx << ", " << std::setw(1) << gy << ") "
                               << std::setw(15) << std::scientific << std::setprecision(6) << real
                               << " "
                               << std::setw(15) << std::scientific << std::setprecision(6) << imag
                               << "\n";
                            ++mode_index;
                        }
                    }
                }

                std::cout << os.str() << std::flush;
            }
        }
    }

    void execute_command(const Command& command) {
        if (command.type == Command::Type::Run) {
            run_steps(command.run.steps);
        }
        else if (command.type == Command::Type::Fix) {
            monitor.print_fix_command(command.fix);
            // TODO: active_fixes.apply(command.fix);
        }
        else if (command.type == Command::Type::Measure) {
            monitor.print_measure_command(*command.measure);
            // TODO: measurements.apply_measure_command(...);
        }
    }

    void run_steps(int nsteps) {
        ++run_index;

        for (int local_step = 1; local_step <= nsteps; ++local_step) {
            ++step;
            monitor.print_progress(run_index, local_step, nsteps, step, time);
        }

        monitor.finish_run_segment(run_index, step, time);
    }

public:
    explicit Project(const Params& p)
        : params(p),
          domain(params),
          state(domain, params),
          buf_physical_state(domain, params),
          fourier(domain),
          registry(build_measure_registry()),
          monitor(params, "output.log", domain.rank() == 0),
          step(0),
          run_index(0),
          time(0.0) {}

    void run() {
        run_fft_sine_mode_test();

        monitor.start();

        for (const auto& command : params.commands) {
            execute_command(command);
        }

        monitor.finish();
    }
};

#endif
