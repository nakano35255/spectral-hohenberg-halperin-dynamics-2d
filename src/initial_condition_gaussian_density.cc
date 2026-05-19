#include "initial_condition_gaussian_density.h"
#include "initial_condition_gaussian_density_style.h"

#include <cmath>
#include <stdexcept>

namespace {
    constexpr double PI = 3.14159265358979323846;

    int signed_ky(int gy, int ny) {
        return gy <= ny / 2 ? gy : gy - ny;
    }
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

void GaussianDensityInitialCondition::apply(State& state, const Domain2D& domain) const {
    const int nx = domain.nx_global();
    const int ny = domain.ny_global();
    const double grid_size = static_cast<double>(nx) * static_cast<double>(ny);
    const double lx = domain.lx();
    const double ly = domain.ly();
    const double gaussian_prefactor =
        amplitude_ * grid_size * 2.0 * PI * sigma_x_ * sigma_y_ / (lx * ly);
    const Box2D& box = domain.spectral_box();

    for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
        const int ky_index = signed_ky(gy, ny);
        const double qy = 2.0 * PI * static_cast<double>(ky_index) / ly;

        for (int kx = box.low[0]; kx <= box.high[0]; ++kx) {
            const double qx = 2.0 * PI * static_cast<double>(kx) / lx;
            const double envelope = std::exp(
                -0.5 * (
                    sigma_x_ * sigma_x_ * qx * qx +
                    sigma_y_ * sigma_y_ * qy * qy
                )
            );
            const double phase = -(qx * x0_ + qy * y0_);

            Complex mode = gaussian_prefactor * envelope * Complex(std::cos(phase), std::sin(phase));

            if (kx == 0 && gy == 0) {
                mode += Complex(base_ * grid_size, 0.0);
            }

            state.rho_hat(component_, kx, gy) = mode;
        }
    }
}
