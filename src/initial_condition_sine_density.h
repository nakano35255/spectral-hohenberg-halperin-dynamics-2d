#ifndef SHHD_INITIAL_CONDITION_SINE_DENSITY_H
#define SHHD_INITIAL_CONDITION_SINE_DENSITY_H

#include "initial_condition.h"

#include <iosfwd>
#include <memory>
#include <string>

class SineDensityInitialCondition : public DensityInitialCondition {
private:
    double base_ = 0.0;
    double amplitude_ = 0.0;
    int nkx_ = 0;
    int nky_ = 0;

public:
    SineDensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command);

    void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const override;
};


#endif