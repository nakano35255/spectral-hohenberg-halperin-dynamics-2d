#include "initial_condition_uniform_density.h"
#include "initial_condition_uniform_density_style.h"


UniformDensityInitialCondition::UniformDensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command) :
    DensityInitialCondition(params, command)
{
    auto cfg = std::dynamic_pointer_cast<const UniformDensityInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("UniformDensityInitialCondition: invalid command type.");
    }

    value_ = cfg->value;
}

void UniformDensityInitialCondition::apply(State& state, const Domain2D& domain) const {
    const double amplitude = value_ * static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    const Box2D& box = domain.spectral_box();

    if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
        state.rho_hat(component_, 0, 0) = Complex(amplitude, 0.0);
    }

}

