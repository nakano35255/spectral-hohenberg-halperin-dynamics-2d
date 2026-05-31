#ifndef SHHD_MODEL_FREE_ENERGY_QUADRATIC_H
#define SHHD_MODEL_FREE_ENERGY_QUADRATIC_H

#include "model_free_energy.h"

#include <vector>

class QuadraticFreeEnergy : public FreeEnergy {
private:
    std::vector<double> coefficients_;

public:
    QuadraticFreeEnergy(
        int num_order_parameters,
        std::vector<double> coefficients
    );

    double chemical_potential_k0_coefficient(int order_parameter) const override;
};

#endif
