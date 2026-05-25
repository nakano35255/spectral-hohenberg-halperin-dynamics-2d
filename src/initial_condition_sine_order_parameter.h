#ifndef SFI_INITIAL_CONDITION_SINE_ORDER_PARAMETER_H
#define SFI_INITIAL_CONDITION_SINE_ORDER_PARAMETER_H

#include "initial_condition.h"

#include <memory>

class SineOrderParameterInitialCondition : public OrderParameterInitialCondition {
private:
    double base_ = 0.0;
    double amplitude_ = 0.0;
    int nkx_ = 0;
    int nky_ = 0;

public:
    SineOrderParameterInitialCondition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif
