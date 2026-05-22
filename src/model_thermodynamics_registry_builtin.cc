#include "model_thermodynamics_registry_builtin.h"
#include "model_thermodynamics_ideal_gas_style.h"


ThermodynamicsRegistry build_thermodynamics_registry() {
    ThermodynamicsRegistry registry;
    registry.register_thermo_style(std::make_unique<IdealGasThermodynamicsStyle>());
    return registry;
}
