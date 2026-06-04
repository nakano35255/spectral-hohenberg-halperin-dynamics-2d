#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_MOMENTUM_INCOMPRESSIBLE_STYLE_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_MOMENTUM_INCOMPRESSIBLE_STYLE_H

#include "initial_condition_equilibrium_gaussian_momentum_incompressible.h"
#include "initial_condition_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace equilibrium_gaussian_incompressible_momentum_initial_condition {
    inline const std::string TYPE_NAME = "equilibrium/gaussian/incompressible";
}

struct EquilibriumGaussianIncompressibleMomentumInitialConditionCommand : public MomentumInitialConditionCommandBase {
    int seed = 12345;
    double kBT = 1.0;

    void print(std::ostream& os) const override {
        const char* component_name = "?";
        if (direction == 0) component_name = "x";
        else if (direction == 1) component_name = "y";

        os << "  " << std::left << std::setw(25)
        << "Momentum" << ": Component " << component_name
        << " " << type
        << " kBT " << kBT
        << " seed " << seed << '\n';
    }
};

class EquilibriumGaussianIncompressibleMomentumInitialConditionStyle : public MomentumInitialConditionStyle {
private:
    const std::string name_ = equilibrium_gaussian_incompressible_momentum_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override { return name_; }

    std::shared_ptr<MomentumInitialConditionCommandBase> parse_command(
        int direction,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        (void)params;

        auto cmd = std::make_shared<EquilibriumGaussianIncompressibleMomentumInitialConditionCommand>();
        cmd->type = name_;
        cmd->direction = direction;
        cmd->seed = std::stoi(args.get_or("seed", "12345"));
        cmd->kBT = std::stod(args.get_required("kBT"));

        if (cmd->kBT < 0.0) {
            throw std::runtime_error("set momentum " + name_ + " requires kBT >= 0.");
        }

        return cmd;
    }

    std::unique_ptr<MomentumInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    ) const override {
        return std::make_unique<EquilibriumGaussianIncompressibleMomentumInitialCondition>(params, command);
    }
};

#endif