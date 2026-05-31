#ifndef SHHD_MODEL_FREE_ENERGY_PHI4_H
#define SHHD_MODEL_FREE_ENERGY_PHI4_H

#include "model_free_energy.h"

#include <vector>

class Phi4FreeEnergy : public FreeEnergy {
private:
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> u_;
    bool has_physical_chemical_potential_ = false;

public:
    Phi4FreeEnergy(
        int num_order_parameters,
        std::vector<double> a,
        std::vector<double> b,
        std::vector<double> u
    );

    double chemical_potential_k0_coefficient(int order_parameter) const override;
    double chemical_potential_k2_coefficient(int order_parameter) const override;

    bool has_physical_chemical_potential() const override;
    double physical_chemical_potential(int order_parameter, const double* psi) const override;

    double physical_free_energy_density(const double* psi) const override;
};

#endif
