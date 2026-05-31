#include "initial_condition_registry_builtin.h"
#include "initial_condition_uniform_density_style.h"
#include "initial_condition_gaussian_density_style.h"
#include "initial_condition_sine_density_style.h"
#include "initial_condition_uniform_momentum_style.h"
#include "initial_condition_sine_momentum_style.h"
#include "initial_condition_random_vorticity_momentum_style.h"
#include "initial_condition_uniform_order_parameter_style.h"
#include "initial_condition_gaussian_order_parameter_style.h"
#include "initial_condition_sine_order_parameter_style.h"


InitialConditionRegistry build_initial_condition_registry() {
    InitialConditionRegistry registry;
    registry.register_density_style(std::make_unique<UniformDensityInitialConditionStyle>());
    registry.register_density_style(std::make_unique<GaussianDensityInitialConditionStyle>());
    registry.register_density_style(std::make_unique<SineDensityInitialConditionStyle>());
    registry.register_momentum_style(std::make_unique<UniformMomentumInitialConditionStyle>());
    registry.register_momentum_style(std::make_unique<SineMomentumInitialConditionStyle>());
    registry.register_momentum_style(std::make_unique<RandomVorticityMomentumInitialConditionStyle>());
    registry.register_order_parameter_style(std::make_unique<UniformOrderParameterInitialConditionStyle>());
    registry.register_order_parameter_style(std::make_unique<GaussianOrderParameterInitialConditionStyle>());
    registry.register_order_parameter_style(std::make_unique<SineOrderParameterInitialConditionStyle>());
    return registry;
}
