#include "model_transport_coefficient_constant.h"

#include <utility>

// ---------------------------------------------------------------------- //
ConstantTransportCoefficient::ConstantTransportCoefficient(
    int num_order_parameters,
    double shear_viscosity,
    double bulk_viscosity,
    std::vector<double> order_parameter_mobility
) : TransportCoefficient(
        num_order_parameters,
        shear_viscosity,
        bulk_viscosity,
        std::move(order_parameter_mobility)
    ) {}
// ---------------------------------------------------------------------- //
