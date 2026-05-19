#ifndef SFI_INITIAL_CONDITION_UNIFORM_DENSITY_H
#define SFI_INITIAL_CONDITION_UNIFORM_DENSITY_H

#include "initial_condition.h"

#include <iosfwd>
#include <memory>
#include <string>

class UniformDensityInitialCondition : public DensityInitialCondition {
private:
    double value_ = 0.0;

public:
    UniformDensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command);

    void apply(State& state, const Domain2D& domain) const override;
};


#endif