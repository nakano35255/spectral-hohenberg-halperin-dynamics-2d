#include "measure_yokota_green_kubo.h"
#include "measure_yokota_green_kubo_style.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    constexpr double PI = 3.14159265358979323846;
}

// ---------------------------------------------------------------------- //
YokotaGreenKuboMeasure::YokotaGreenKuboMeasure(
    const Params& params,
    const Domain2D& domain,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    spectral_mask_(params, domain)
{
    auto cfg = std::dynamic_pointer_cast<const YokotaGreenKuboMeasureCommand>(command);
    if (!cfg) {
        throw std::runtime_error("YokotaGreenKuboMeasure: invalid command type.");
    }

    nevery_ = cfg->nevery;
    nblock_ = cfg->nblock;
    file_ = cfg->file;
    mode_selection_ = cfg->mode;
    kBT_ = cfg->kBT;

    configure_modes();
}
// ---------------------------------------------------------------------- //
FluxRequest YokotaGreenKuboMeasure::flux_request() const {
    FluxRequest request;
    request.momentum = true;
    return request;
}
// ---------------------------------------------------------------------- //
int YokotaGreenKuboMeasure::wrap_gy(int gy) const {
    const int ny = domain_.ny_global();
    int wrapped = gy % ny;
    if (wrapped < 0) {
        wrapped += ny;
    }
    return wrapped;
}
// ---------------------------------------------------------------------- //
double YokotaGreenKuboMeasure::cell_area() const {
    return domain_.lx() * domain_.ly()
        / (static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global()));
}
// ---------------------------------------------------------------------- //
double YokotaGreenKuboMeasure::system_area() const {
    return domain_.lx() * domain_.ly();
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::add_mode(int requested_nx, int requested_ny) {
    const int nx = domain_.nx_global();
    const int ny = domain_.ny_global();

    if (std::abs(requested_nx) > nx / 2) {
        throw std::runtime_error("YokotaGreenKuboMeasure: requested nx is outside the spectral range.");
    }
    if (std::abs(requested_ny) > ny / 2) {
        throw std::runtime_error("YokotaGreenKuboMeasure: requested ny is outside the spectral range.");
    }
    if (requested_nx == 0 && requested_ny == 0) {
        return;
    }

    ModeAccumulator mode;
    mode.requested_nx = requested_nx;
    mode.requested_ny = requested_ny;

    if (requested_nx < 0) {
        mode.stored_gx = -requested_nx;
        mode.stored_gy = wrap_gy(-requested_ny);
        mode.use_conjugate = true;
    }
    else {
        mode.stored_gx = requested_nx;
        mode.stored_gy = wrap_gy(requested_ny);
        mode.use_conjugate = false;
    }

    if (!spectral_mask_.active(mode.stored_gx, mode.stored_gy)) {
        throw std::runtime_error("YokotaGreenKuboMeasure: requested wave vector is outside the active spectral mask.");
    }

    mode.kx = 2.0 * PI * static_cast<double>(requested_nx) / domain_.lx();
    mode.ky = 2.0 * PI * static_cast<double>(requested_ny) / domain_.ly();

    const Box2D& box = domain_.spectral_box();
    mode.owns_mode = !box.empty()
                  && mode.stored_gx >= box.low[0] && mode.stored_gx <= box.high[0]
                  && mode.stored_gy >= box.low[1] && mode.stored_gy <= box.high[1];

    if (mode.owns_mode) {
        mode.mode_index = domain_.spectral_local_index(mode.stored_gx, mode.stored_gy);
    }

    const int ntau = nblock_ / nevery_;
    mode.sum_by_tau.assign(static_cast<std::size_t>(ntau), 0.0);
    mode.count_by_tau.assign(static_cast<std::size_t>(ntau), 0);

    modes_.push_back(std::move(mode));
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::configure_modes() {
    const int max_nx = params_.grid.active_num_grid[0] / 2 - 1;
    const int max_ny = params_.grid.active_num_grid[1] / 2 - 1;

    if (mode_selection_ == "diagonal") {
        const int max_n = std::min(max_nx, max_ny);
        for (int n = 1; n <= max_n; ++n) {
            add_mode(n, n);
        }
    }
    else if (mode_selection_ == "all") {
        for (int nx = 0; nx <= max_nx; ++nx) {
            for (int ny = -max_ny; ny <= max_ny; ++ny) {
                if (nx == 0 && ny <= 0) {
                    continue;
                }
                add_mode(nx, ny);
            }
        }
    }
    else {
        throw std::runtime_error("YokotaGreenKuboMeasure: mode must be diagonal|all.");
    }

    if (modes_.empty()) {
        throw std::runtime_error("YokotaGreenKuboMeasure: no active wave vector is selected.");
    }
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::open_output_if_needed() {
    int open_ok = 1;

    if (domain_.rank() == 0 && !output_.is_open()) {
        output_.open(file_);
        if (!output_) {
            open_ok = 0;
        }
        else {
            output_ << "# step time block tau mode_index nx ny kx ky S_block S_mean nsample integral_real integral_imag\n";
            output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain_.comm());
    if (open_ok == 0) {
        throw std::runtime_error("YokotaGreenKuboMeasure: cannot open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::observe(
    const State& /*state*/,
    FourierTransform2D& /*fft*/,
    MeasureWorkspace& /*workspace*/,
    const FluxBuffer& flux,
    int step,
    double time
) {
    open_output_if_needed();

    ++block_step_;

    const Complex* pi_xy = flux.momentum_flux_xy_hat_data();
    const double endpoint_weight = (block_step_ == 1) ? 0.5 : 1.0;
    const double integral_weight = endpoint_weight * params_.runtime.dt * cell_area();

    for (ModeAccumulator& mode : modes_) {
        if (!mode.owns_mode) {
            continue;
        }

        Complex value = pi_xy[mode.mode_index];
        if (mode.use_conjugate) {
            value = std::conj(value);
        }

        mode.integral += integral_weight * value;
    }

    if (block_step_ % nevery_ == 0) {
        const double tau = static_cast<double>(block_step_) * params_.runtime.dt;
        const int tau_index = block_step_ / nevery_ - 1;
        const std::size_t nmodes = modes_.size();

        std::vector<double> local_integrals(2 * nmodes, 0.0);
        for (std::size_t i = 0; i < nmodes; ++i) {
            local_integrals[2 * i] = modes_[i].integral.real();
            local_integrals[2 * i + 1] = modes_[i].integral.imag();
        }

        std::vector<double> global_integrals(2 * nmodes, 0.0);
        MPI_Reduce(
            local_integrals.data(), global_integrals.data(),
            static_cast<int>(global_integrals.size()), MPI_DOUBLE, MPI_SUM,
            0, domain_.comm()
        );

        int write_ok = 1;

        if (domain_.rank() == 0) {
            for (std::size_t i = 0; i < nmodes; ++i) {
                ModeAccumulator& mode = modes_[i];
                const Complex total_integral(global_integrals[2 * i], global_integrals[2 * i + 1]);
                const double s_block = std::norm(total_integral) / (2.0 * tau * kBT_ * system_area());

                const std::size_t index = static_cast<std::size_t>(tau_index);
                mode.sum_by_tau[index] += s_block;
                ++mode.count_by_tau[index];

                const double s_mean =
                    mode.sum_by_tau[index] / static_cast<double>(mode.count_by_tau[index]);

                output_ << step
                        << ' ' << time
                        << ' ' << block_index_
                        << ' ' << tau
                        << ' ' << i
                        << ' ' << mode.requested_nx
                        << ' ' << mode.requested_ny
                        << ' ' << mode.kx
                        << ' ' << mode.ky
                        << ' ' << s_block
                        << ' ' << s_mean
                        << ' ' << mode.count_by_tau[index]
                        << ' ' << total_integral.real()
                        << ' ' << total_integral.imag()
                        << '\n';
            }

            if (!output_) {
                write_ok = 0;
            }
        }

        MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
        if (write_ok == 0) {
            throw std::runtime_error("YokotaGreenKuboMeasure: failed to write file: " + file_);
        }
    }

    if (block_step_ == nblock_) {
        ++block_index_;
        block_step_ = 0;
        for (ModeAccumulator& mode : modes_) {
            mode.integral = Complex(0.0, 0.0);
        }
    }
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::finalize() {
    if (output_.is_open()) {
        output_.close();
    }
}
// ---------------------------------------------------------------------- //
