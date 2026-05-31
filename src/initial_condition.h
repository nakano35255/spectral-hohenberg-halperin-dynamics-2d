#ifndef SHHD_INITIAL_CONDITION_H
#define SHHD_INITIAL_CONDITION_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "spectral_mask.h"

class DensityInitialCondition {
private:
    std::shared_ptr<const DensityInitialConditionCommandBase> command_;

public:
    DensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command)
        : command_(std::move(command))
    {
        (void)params;
    }

    virtual ~DensityInitialCondition() = default;

    virtual void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const = 0;
};

class MomentumInitialCondition {
private:
    std::shared_ptr<const MomentumInitialConditionCommandBase> command_;

protected:
    int direction_ = -1;

public:
    MomentumInitialCondition(const Params& params, std::shared_ptr<const MomentumInitialConditionCommandBase> command)
        : command_(std::move(command)),
          direction_(command_->direction)
    {
        (void)params;
    }

    virtual ~MomentumInitialCondition() = default;

    virtual void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const = 0;

};

class OrderParameterInitialCondition {
private:
    std::shared_ptr<const OrderParameterInitialConditionCommandBase> command_;

protected:
    int component_ = -1;

public:
    OrderParameterInitialCondition(const Params& params, std::shared_ptr<const OrderParameterInitialConditionCommandBase> command)
        : command_(std::move(command)),
          component_(command_->component)
    {
        (void)params;
    }

    virtual ~OrderParameterInitialCondition() = default;

    virtual void apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const = 0;

};

void validate_initial_momentum_is_transverse(const State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask);

#endif
