#include "measure_yokota_green_kubo.h"
#include "measure_yokota_green_kubo_style.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
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
    kBT_ = params.fix.noise.kBT;

    if (kBT_ <= 0.0) {
        throw std::runtime_error("YokotaGreenKuboMeasure requires positive kBT from fix noise.");
    }

    max_nx_ = std::max(0, params_.grid.active_num_grid[0] / 2 - 1);
    max_ny_ = std::max(0, params_.grid.active_num_grid[1] / 2 - 1);
    ntau_ = nblock_ / nevery_;

    if (mode_selection_ == "diagonal") {
        nmodes_ = static_cast<std::size_t>(std::min(max_nx_, max_ny_));
    }
    else if (mode_selection_ == "all") {
        const int ny_count = 2 * max_ny_ + 1;
        nmodes_ = static_cast<std::size_t>(max_ny_)
                + static_cast<std::size_t>(max_nx_) * static_cast<std::size_t>(ny_count);
    }
    else {
        throw std::runtime_error("YokotaGreenKuboMeasure: mode must be diagonal|all.");
    }

    if (nmodes_ == 0) {
        throw std::runtime_error("YokotaGreenKuboMeasure: no active wave vector is selected.");
    }

    owned_local_indices_.assign(nmodes_, std::nullopt);
    integrals_.assign(nmodes_, Complex(0.0, 0.0));
    sum_by_mode_tau_.assign(nmodes_ * static_cast<std::size_t>(ntau_), 0.0);

    for (std::size_t mode_index = 0; mode_index < nmodes_; ++mode_index) {
        const auto [nx, ny] = wave_number(mode_index);
        owned_local_indices_[mode_index] = local_mode_index(nx, ny);
    }
}
// ---------------------------------------------------------------------- //
std::pair<int, int> YokotaGreenKuboMeasure::wave_number(std::size_t mode_index) const {
    if (mode_selection_ == "diagonal") {
        const int n = static_cast<int>(mode_index) + 1;
        return {n, n};
    }

    const std::size_t y_axis_count = static_cast<std::size_t>(max_ny_);
    if (mode_index < y_axis_count) {
        return {0, static_cast<int>(mode_index) + 1};
    }

    const std::size_t shifted = mode_index - y_axis_count;
    const int ny_count = 2 * max_ny_ + 1;

    const int nx = static_cast<int>(shifted / static_cast<std::size_t>(ny_count)) + 1;
    const int ny = static_cast<int>(shifted % static_cast<std::size_t>(ny_count)) - max_ny_;

    return {nx, ny};
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
std::optional<std::size_t> YokotaGreenKuboMeasure::local_mode_index(int nx, int ny) const {
    const int gx = nx;
    int gy = wrap_gy(ny);

    if (!spectral_mask_.active(gx, gy)) {
        throw std::runtime_error("YokotaGreenKuboMeasure: selected wave vector is outside the active spectral mask.");
    }

    const Box2D& box = domain_.spectral_box();
    const bool owns_mode = !box.empty()
                        && gx >= box.low[0] && gx <= box.high[0] && gy >= box.low[1] && gy <= box.high[1];

    if (!owns_mode) {
        return std::nullopt;
    }

    return domain_.spectral_local_index(gx, gy);
}
// ---------------------------------------------------------------------- //
std::size_t YokotaGreenKuboMeasure::mode_tau_index(std::size_t mode_index, int tau_index) const {
    return mode_index * static_cast<std::size_t>(ntau_) + static_cast<std::size_t>(tau_index);
}
// ---------------------------------------------------------------------- //

FluxRequest YokotaGreenKuboMeasure::flux_request() const {
    FluxRequest request;
    request.momentum = true;
    return request;
}
// ---------------------------------------------------------------------- //
void YokotaGreenKuboMeasure::observe(
    const State& /*state*/,
    FourierTransform2D& /*fft*/,
    MeasureWorkspace& /*workspace*/,
    const FluxBuffer& flux,
    int /*step*/,
    double /*time*/
) {
    ++block_step_;

    const double area = domain_.lx() * domain_.ly();
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());
    const double cell_area = area / grid_size;

    const Complex* pi_xy = flux.momentum_flux_xy_hat_data();
    const double weighted_dt = (block_step_ == 1) ? 0.5 * params_.runtime.dt : params_.runtime.dt;
    const double integral_weight = weighted_dt * cell_area;

    for (std::size_t mode_index = 0; mode_index < nmodes_; ++mode_index) {
        const auto local_index = owned_local_indices_[mode_index];
        if (!local_index) {
            continue;
        }

        integrals_[mode_index] += integral_weight * pi_xy[*local_index];
    }

    if (block_step_ % nevery_ == 0) {
        const int tau_index = block_step_ / nevery_ - 1;
        const double tau = static_cast<double>(block_step_) * params_.runtime.dt;

        std::vector<double> local_integrals(2 * nmodes_, 0.0);
        for (std::size_t mode_index = 0; mode_index < nmodes_; ++mode_index) {
            local_integrals[2 * mode_index] = integrals_[mode_index].real();
            local_integrals[2 * mode_index + 1] = integrals_[mode_index].imag();
        }

        std::vector<double> global_integrals(2 * nmodes_, 0.0);
        MPI_Reduce(local_integrals.data(), global_integrals.data(), static_cast<int>(global_integrals.size()), MPI_DOUBLE, MPI_SUM, 0, domain_.comm());

        if (domain_.rank() == 0) {
            for (std::size_t mode_index = 0; mode_index < nmodes_; ++mode_index) {
                const Complex total_integral(global_integrals[2 * mode_index], global_integrals[2 * mode_index + 1]);
                const double s_block = std::norm(total_integral) / (2.0 * tau * kBT_ * area);
                sum_by_mode_tau_[mode_tau_index(mode_index, tau_index)] += s_block;
            }
        }
    }

    if (block_step_ == nblock_) {
        ++nsample_;

        int write_ok = 1;
        if (domain_.rank() == 0) {
            write_ok = write_output_table() ? 1 : 0;
        }

        MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
        if (write_ok == 0) {
            throw std::runtime_error("YokotaGreenKuboMeasure: failed to write file: " + file_);
        }

        block_step_ = 0;
        std::fill(integrals_.begin(), integrals_.end(), Complex(0.0, 0.0));
    }
}
// ---------------------------------------------------------------------- //
bool YokotaGreenKuboMeasure::write_output_table() const {
    std::ofstream out(file_, std::ios::out | std::ios::trunc);
    if (!out) {
        return false;
    }

    out << "# nsample tau mode_index kx ky S_mean\n";
    out << std::scientific << std::setprecision(16);

    for (int tau_index = 0; tau_index < ntau_; ++tau_index) {
        const double tau =
            static_cast<double>((tau_index + 1) * nevery_) * params_.runtime.dt;

        for (std::size_t mode_index = 0; mode_index < nmodes_; ++mode_index) {
            const auto [nx, ny] = wave_number(mode_index);
            const double kx = 2.0 * PI * static_cast<double>(nx) / domain_.lx();
            const double ky = 2.0 * PI * static_cast<double>(ny) / domain_.ly();
            const double s_mean = sum_by_mode_tau_[mode_tau_index(mode_index, tau_index)] / static_cast<double>(nsample_);

            out << nsample_
                << ' ' << tau
                << ' ' << mode_index
                << ' ' << kx
                << ' ' << ky
                << ' ' << s_mean
                << '\n';
        }
    }

    return static_cast<bool>(out);
}
