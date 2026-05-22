#include "model_transport_coefficient_registry_builtin.h"
#include "model_transport_coefficient_constant_style.h"


TransportCoefficientRegistry build_transport_coefficient_registry() {
    TransportCoefficientRegistry registry;
    registry.register_transport_style(std::make_unique<ConstantTransportCoefficientStyle>());
    return registry;
}
