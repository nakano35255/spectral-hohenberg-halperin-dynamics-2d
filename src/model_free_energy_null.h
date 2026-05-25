#ifndef SFI_MODEL_FREE_ENERGY_NULL_H
#define SFI_MODEL_FREE_ENERGY_NULL_H

#include "model_free_energy.h"

class NullFreeEnergy : public FreeEnergy {
public:
    explicit NullFreeEnergy(int num_order_parameters);
};

#endif
