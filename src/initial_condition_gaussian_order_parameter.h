#ifndef SHHD_INITIAL_CONDITION_GAUSSIAN_ORDER_PARAMETER_H
#define SHHD_INITIAL_CONDITION_GAUSSIAN_ORDER_PARAMETER_H

#include "initial_condition.h"

#include <memory>

class GaussianOrderParameterInitialCondition : public OrderParameterInitialCondition {
private:
    double base_ = 0.0;
    double amplitude_ = 0.0;
    double x0_ = 0.0;
    double y0_ = 0.0;
    double sigma_x_ = 1.0;
    double sigma_y_ = 1.0;

public:
    GaussianOrderParameterInitialCondition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif
