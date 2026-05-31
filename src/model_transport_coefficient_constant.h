#ifndef SHHD_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_H
#define SHHD_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_H

#include "model_transport_coefficient.h"

#include <vector>

class ConstantTransportCoefficient : public TransportCoefficient {
public:
    explicit ConstantTransportCoefficient(
        int num_order_parameters,
        double shear_viscosity,
        double bulk_viscosity,
        std::vector<double> order_parameter_mobility
    );
};

#endif
