#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_DENSITY_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_DENSITY_H

#include "initial_condition.h"

class EquilibriumGaussianDensityInitialCondition : public DensityInitialCondition {
private:
    int seed_ = 12345;
    double kBT_ = 1.0;
    double rho0_ = 1.0;
    double pressure_coefficient_ = 0.0;

public:
    EquilibriumGaussianDensityInitialCondition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif
