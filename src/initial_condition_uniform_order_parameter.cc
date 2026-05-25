#include "initial_condition_uniform_order_parameter.h"
#include "initial_condition_uniform_order_parameter_style.h"

#include <stdexcept>

UniformOrderParameterInitialCondition::UniformOrderParameterInitialCondition(
    const Params& params,
    std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
)
    : OrderParameterInitialCondition(params, command)
{
    const auto cfg = std::dynamic_pointer_cast<const UniformOrderParameterInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("UniformOrderParameterInitialCondition: invalid command type.");
    }

    value_ = cfg->value;
}

void UniformOrderParameterInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    const double amplitude = value_ * static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    Complex* psi = state.psi_hat_data(component_);

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            psi[mode.index] = Complex(amplitude, 0.0);
            break;
        }
    }
}
