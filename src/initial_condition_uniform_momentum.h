#ifndef SFI_INITIAL_CONDITION_UNIFORM_MOMENTUM_H
#define SFI_INITIAL_CONDITION_UNIFORM_MOMENTUM_H

#include "initial_condition.h"

#include <iosfwd>
#include <memory>
#include <string>

class UniformMomentumInitialCondition : public MomentumInitialCondition {
private:
    double value_ = 0.0;

public:
    UniformMomentumInitialCondition(const Params& params, std::shared_ptr<const MomentumInitialConditionCommandBase> command);

    void apply(State& state, const Domain2D& domain) const override;
};


#endif