#include "model_free_energy_registry_builtin.h"
#include "model_free_energy_phi4_style.h"
#include "model_free_energy_quadratic_style.h"

FreeEnergyRegistry build_free_energy_registry() {
    FreeEnergyRegistry registry;
    registry.register_free_energy_style(std::make_unique<QuadraticFreeEnergyStyle>());
    registry.register_free_energy_style(std::make_unique<Phi4FreeEnergyStyle>());
    return registry;
}
