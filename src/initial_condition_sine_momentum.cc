#include "initial_condition_sine_momentum.h"
#include "initial_condition_sine_momentum_style.h"

#include <stdexcept>

namespace {
    int wrap_ky(int nky, int ny) {
        return nky >= 0 ? nky : ny + nky;
    }
}

SineMomentumInitialCondition::SineMomentumInitialCondition(
    const Params& params,
    std::shared_ptr<const MomentumInitialConditionCommandBase> command
) : MomentumInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const SineMomentumInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("SineMomentumInitialCondition: invalid command type.");
    }

    base_ = cfg->base;
    amplitude_ = cfg->amplitude;
    nkx_ = cfg->nkx;
    nky_ = cfg->nky;
}

void SineMomentumInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    const int nx = domain.nx_global();
    const int ny = domain.ny_global();
    const double grid_size = static_cast<double>(nx) * static_cast<double>(ny);

    const int ky = wrap_ky(nky_, ny);
    const int conjugate_ky = (ny - ky) % ny;
    const Complex positive_mode(0.0, -0.5 * amplitude_ * grid_size);

    Complex* momentum = nullptr;
    if (direction_ == 0) momentum = state.jx_hat_data();
    else if (direction_ == 1) momentum = state.jy_hat_data();
    else throw std::runtime_error("SineMomentumInitialCondition: direction must be 0 or 1.");

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            momentum[mode.index] = Complex(base_ * grid_size, 0.0);
        }

        if (mode.gx == nkx_ && mode.gy == ky) {
            momentum[mode.index] += positive_mode;
        }

        if (conjugate_ky != ky && mode.gx == nkx_ && mode.gy == conjugate_ky) {
            momentum[mode.index] += std::conj(positive_mode);
        }
    }
}
