#ifndef SFI_INITIAL_CONDITION_UNIFORM_ORDER_PARAMETER_H
#define SFI_INITIAL_CONDITION_UNIFORM_ORDER_PARAMETER_H

#include "initial_condition.h"

#include <memory>

class UniformOrderParameterInitialCondition : public OrderParameterInitialCondition {
private:
    double value_ = 0.0;

public:
    UniformOrderParameterInitialCondition(const Params& params, std::shared_ptr<const OrderParameterInitialConditionCommandBase> command);

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif