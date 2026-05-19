#include "initial_condition_registry_builtin.h"
#include "initial_condition_uniform_density_style.h"
#include "initial_condition_uniform_momentum_style.h"

InitialConditionRegistry build_initial_condition_registry() {
    InitialConditionRegistry registry;
    registry.register_density_style(std::make_unique<UniformDensityInitialConditionStyle>());
    registry.register_momentum_style(std::make_unique<UniformMomentumInitialConditionStyle>());
    return registry;
}
