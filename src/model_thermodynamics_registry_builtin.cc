#include "model_thermodynamics_registry_builtin.h"
#include "model_thermodynamics_ideal_gas_style.h"


ThermodynamicsModelRegistry build_thermodynamics_model_registry() {
    ThermodynamicsModelRegistry registry;
    registry.register_thermo_style(std::make_unique<IdealGasThermodynamicsModelStyle>());
    return registry;
}
