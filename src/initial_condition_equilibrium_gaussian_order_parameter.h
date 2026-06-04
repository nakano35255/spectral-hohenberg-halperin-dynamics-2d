#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_ORDER_PARAMETER_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_ORDER_PARAMETER_H

#include "initial_condition.h"

class EquilibriumGaussianOrderParameterInitialCondition : public OrderParameterInitialCondition {
private:
    int seed_ = 12345;
    double kBT_ = 1.0;
    double psi0_ = 0.0;
    double k0_coefficient_ = 0.0;
    double k2_coefficient_ = 0.0;
    double k4_coefficient_ = 0.0;

public:
    EquilibriumGaussianOrderParameterInitialCondition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};

#endif
