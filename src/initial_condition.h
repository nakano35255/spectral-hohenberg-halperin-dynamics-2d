#ifndef SFI_INITIAL_CONDITION_H
#define SFI_INITIAL_CONDITION_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"

class DensityInitialCondition {
private:
    std::shared_ptr<const DensityInitialConditionCommandBase> command_;

protected:
    int component_ = -1;

public:
    DensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command)
        : command_(std::move(command)),
          component_(command_->component) {
        (void)params;
    }

    virtual ~DensityInitialCondition() = default;

    virtual void apply(State& state, const Domain2D& domain) const = 0;
};

class MomentumInitialCondition {
private:
    std::shared_ptr<const MomentumInitialConditionCommandBase> command_;

protected:
    int component_ = -1;

public:
     MomentumInitialCondition(const Params& params, std::shared_ptr<const MomentumInitialConditionCommandBase> command)
        : command_(std::move(command)),
          component_(command_->component) {
        (void)params;
    }

    virtual ~MomentumInitialCondition() = default;

    virtual void apply(State& state, const Domain2D& domain) const = 0;

};

#endif
