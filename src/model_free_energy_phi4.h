#ifndef SFI_MODEL_FREE_ENERGY_PHI4_H
#define SFI_MODEL_FREE_ENERGY_PHI4_H

#include "model_free_energy.h"

#include <vector>

class Phi4FreeEnergy : public FreeEnergy {
private:
    std::vector<double> k0_;
    std::vector<double> k2_;
    std::vector<double> phi4_;
    bool has_physical_chemical_potential_ = false;

public:
    Phi4FreeEnergy(
        int num_order_parameters,
        std::vector<double> k0,
        std::vector<double> k2,
        std::vector<double> phi4
    );

    double chemical_potential_k0_coefficient(int order_parameter) const override;
    double chemical_potential_k2_coefficient(int order_parameter) const override;

    bool has_physical_chemical_potential() const override;
    double physical_chemical_potential(int order_parameter, const double* psi) const override;
};

#endif
