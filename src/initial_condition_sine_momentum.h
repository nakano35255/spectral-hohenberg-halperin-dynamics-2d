#ifndef SFI_INITIAL_CONDITION_SINE_MOMENTUM_H
#define SFI_INITIAL_CONDITION_SINE_MOMENTUM_H

#include "initial_condition.h"

#include <memory>

class SineMomentumInitialCondition : public MomentumInitialCondition {
private:
    double base_ = 0.0;
    double amplitude_ = 0.0;
    int nkx_ = 0;
    int nky_ = 0;

public:
    SineMomentumInitialCondition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain) const override;
};

#endif
