#include "initial_condition_sine_momentum.h"
#include "initial_condition_sine_momentum_style.h"

#include <stdexcept>

namespace {
    bool owns_mode(const Box2D& box, int kx, int ky) {
        return box.low[0] <= kx && kx <= box.high[0] &&
               box.low[1] <= ky && ky <= box.high[1];
    }

    int wrap_ky(int nky, int ny) {
        return nky >= 0 ? nky : ny + nky;
    }

    Complex& momentum_hat(State& state, int component, int kx, int ky) {
        if (component == 0) {
            return state.jx_hat(kx, ky);
        }
        if (component == 1) {
            return state.jy_hat(kx, ky);
        }
        throw std::runtime_error("SineMomentumInitialCondition: component must be 0 or 1.");
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

void SineMomentumInitialCondition::apply(State& state, const Domain2D& domain) const {
    const int nx = domain.nx_global();
    const int ny = domain.ny_global();
    const double grid_size = static_cast<double>(nx) * static_cast<double>(ny);
    const Box2D& box = domain.spectral_box();

    if (owns_mode(box, 0, 0)) {
        momentum_hat(state, component_, 0, 0) = Complex(base_ * grid_size, 0.0);
    }

    const int ky = wrap_ky(nky_, ny);
    const Complex positive_mode(0.0, -0.5 * amplitude_ * grid_size);

    if (owns_mode(box, nkx_, ky)) {
        momentum_hat(state, component_, nkx_, ky) += positive_mode;
    }

    if (nkx_ == 0 || nkx_ == nx / 2) {
        const int conjugate_ky = (ny - ky) % ny;
        if (conjugate_ky != ky) {
            if (owns_mode(box, nkx_, conjugate_ky)) {
                momentum_hat(state, component_, nkx_, conjugate_ky) += std::conj(positive_mode);
            }
        }
    }
}
