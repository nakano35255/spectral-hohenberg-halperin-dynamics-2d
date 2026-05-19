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

void UniformMomentumInitialCondition::apply(State& state, const Domain2D& domain) const {
    const double amplitude = value_ * static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    const Box2D& box = domain.spectral_box();

    if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
        if (component_ == 0) {
            state.jx_hat(0, 0) = Complex(amplitude, 0.0);
        } else if (component_ == 1) {
            state.jy_hat(0, 0) = Complex(amplitude, 0.0);
        } else {
            throw std::runtime_error("UniformMomentumInitialCondition: component must be 0 or 1.");
        }
    }

}
