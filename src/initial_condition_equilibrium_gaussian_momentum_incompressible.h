#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_MOMENTUM_INCOMPRESSIBLE_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_MOMENTUM_INCOMPRESSIBLE_H

#include "initial_condition.h"

class EquilibriumGaussianIncompressibleMomentumInitialCondition : public MomentumInitialCondition {
private:
    int seed_ = 12345;
    double kBT_ = 1.0;

public:
    EquilibriumGaussianIncompressibleMomentumInitialCondition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif