#ifndef SHHD_INITIAL_CONDITION_RANDOM_VORTICITY_MOMENTUM_H
#define SHHD_INITIAL_CONDITION_RANDOM_VORTICITY_MOMENTUM_H

#include "initial_condition.h"

#include <memory>

class RandomVorticityMomentumInitialCondition : public MomentumInitialCondition {
private:
    double n0_ = 0.0;
    double sigma_ = 0.0;
    double urms_ = 0.0;
    int seed_ = 0;

public:
    RandomVorticityMomentumInitialCondition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif
