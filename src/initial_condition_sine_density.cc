#include "initial_condition_sine_density.h"
#include "initial_condition_sine_density_style.h"

#include <stdexcept>

namespace {
    bool owns_mode(const Box2D& box, int kx, int ky) {
        return box.low[0] <= kx && kx <= box.high[0] &&
               box.low[1] <= ky && ky <= box.high[1];
    }

    int wrap_ky(int nky, int ny) {
        return nky >= 0 ? nky : ny + nky;
    }
}

SineDensityInitialCondition::SineDensityInitialCondition(
    const Params& params,
    std::shared_ptr<const DensityInitialConditionCommandBase> command
) : DensityInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const SineDensityInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("SineDensityInitialCondition: invalid command type.");
    }

    base_ = cfg->base;
    amplitude_ = cfg->amplitude;
    nkx_ = cfg->nkx;
    nky_ = cfg->nky;
}

void SineDensityInitialCondition::apply(State& state, const Domain2D& domain) const {
    const int nx = domain.nx_global();
    const int ny = domain.ny_global();
    const double grid_size = static_cast<double>(nx) * static_cast<double>(ny);
    const Box2D& box = domain.spectral_box();

    if (owns_mode(box, 0, 0)) {
        state.rho_hat(component_, 0, 0) = Complex(base_ * grid_size, 0.0);
    }

    const int ky = wrap_ky(nky_, ny);
    const Complex positive_mode(0.0, -0.5 * amplitude_ * grid_size);

    if (owns_mode(box, nkx_, ky)) {
        state.rho_hat(component_, nkx_, ky) += positive_mode;
    }

    if (nkx_ == 0 || nkx_ == nx / 2) {
        const int conjugate_ky = (ny - ky) % ny;
        if (conjugate_ky != ky) {
            if (owns_mode(box, nkx_, conjugate_ky)) {
                state.rho_hat(component_, nkx_, conjugate_ky) += std::conj(positive_mode);
            }
        }
    }
}
