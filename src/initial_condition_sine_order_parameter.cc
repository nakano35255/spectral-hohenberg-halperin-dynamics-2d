#include "initial_condition_sine_order_parameter.h"
#include "initial_condition_sine_order_parameter_style.h"

#include <algorithm>
#include <stdexcept>

namespace {
    int wrap_ky(int nky, int ny) {
        return nky >= 0 ? nky : ny + nky;
    }
}

SineOrderParameterInitialCondition::SineOrderParameterInitialCondition(
    const Params& params,
    std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
) : OrderParameterInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const SineOrderParameterInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("SineOrderParameterInitialCondition: invalid command type.");
    }

    base_ = cfg->base;
    amplitude_ = cfg->amplitude;
    nkx_ = cfg->nkx;
    nky_ = cfg->nky;
}

void SineOrderParameterInitialCondition::apply(
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

    Complex* psi = state.psi_hat_data(component_);
    std::fill(psi, psi + domain.spectral_size(), Complex(0.0, 0.0));

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            psi[mode.index] = Complex(base_ * grid_size, 0.0);
        }

        if (mode.gx == nkx_ && mode.gy == ky) {
            psi[mode.index] = positive_mode;
        }

        if (conjugate_ky != ky && mode.gx == nkx_ && mode.gy == conjugate_ky) {
            psi[mode.index] = std::conj(positive_mode);
        }
    }
}
