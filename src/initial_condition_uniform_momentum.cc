#include "initial_condition_uniform_momentum.h"
#include "initial_condition_uniform_momentum_style.h"

#include <stdexcept>


UniformMomentumInitialCondition::UniformMomentumInitialCondition(const Params& params, std::shared_ptr<const MomentumInitialConditionCommandBase> command) :
    MomentumInitialCondition(params, command)
{
    auto cfg = std::dynamic_pointer_cast<const UniformMomentumInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("UniformMomentumInitialCondition: invalid command type.");
    }

    value_ = cfg->value;
}

void UniformMomentumInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    const double amplitude = value_ * static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    Complex* momentum = nullptr;
    if (direction_ == 0) momentum = state.jx_hat_data();
    else if (direction_ == 1) momentum = state.jy_hat_data();
    else throw std::runtime_error("UniformMomentumInitialCondition: direction must be 0 or 1.");

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            momentum[mode.index] = Complex(amplitude, 0.0);
            break;
        }
    }

}
