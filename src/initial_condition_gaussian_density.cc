#include "initial_condition_gaussian_density.h"
#include "initial_condition_gaussian_density_style.h"

#include <cmath>
#include <stdexcept>

namespace {
    constexpr double PI = 3.14159265358979323846;
}

GaussianDensityInitialCondition::GaussianDensityInitialCondition(
    const Params& params,
    std::shared_ptr<const DensityInitialConditionCommandBase> command
) : DensityInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const GaussianDensityInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("GaussianDensityInitialCondition: invalid command type.");
    }

    base_ = cfg->base;
    amplitude_ = cfg->amplitude;
    x0_ = cfg->x0;
    y0_ = cfg->y0;
    sigma_x_ = cfg->sigma_x;
    sigma_y_ = cfg->sigma_y;
}

void GaussianDensityInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    const double grid_size =
        static_cast<double>(domain.nx_global()) *
        static_cast<double>(domain.ny_global());

    const double lx = domain.lx();
    const double ly = domain.ly();

    const double gaussian_prefactor = amplitude_ * grid_size * 2.0 * PI * sigma_x_ * sigma_y_ / (lx * ly);

    Complex* rho = state.rho_hat_data();

    for (const SpectralMode2D& spectral_mode : spectral_mask.active_modes()) {
        const int gx = spectral_mode.gx;
        const int gy = spectral_mode.gy;

        const double kx = spectral_mode.kx;
        const double ky = spectral_mode.ky;

        const double envelope = std::exp(
            -0.5 * (sigma_x_ * sigma_x_ * kx * kx + sigma_y_ * sigma_y_ * ky * ky)
        );

        const double phase = -(kx * x0_ + ky * y0_);

        Complex mode = gaussian_prefactor * envelope * Complex(std::cos(phase), std::sin(phase));

        if (gx == 0 && gy == 0) {
            mode += Complex(base_ * grid_size, 0.0);
        }

        rho[spectral_mode.index] = mode;
    }
}
