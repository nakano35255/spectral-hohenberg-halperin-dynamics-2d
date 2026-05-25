#include "model_thermodynamics_registry_builtin.h"
#include "model_thermodynamics_linear_eos_style.h"


ThermodynamicsRegistry build_thermodynamics_registry() {
    ThermodynamicsRegistry registry;
    registry.register_thermo_style(std::make_unique<LinearEOSThermodynamicsStyle>());
    return registry;
}
