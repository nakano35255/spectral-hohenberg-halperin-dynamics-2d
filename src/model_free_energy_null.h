#ifndef SHHD_MODEL_FREE_ENERGY_NULL_H
#define SHHD_MODEL_FREE_ENERGY_NULL_H

#include "model_free_energy.h"

class NullFreeEnergy : public FreeEnergy {
public:
    explicit NullFreeEnergy(int num_order_parameters);
};

#endif
