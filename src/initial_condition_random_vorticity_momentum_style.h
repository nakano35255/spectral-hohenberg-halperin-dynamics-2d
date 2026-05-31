#ifndef SHHD_INITIAL_CONDITION_RANDOM_VORTICITY_MOMENTUM_STYLE_H
#define SHHD_INITIAL_CONDITION_RANDOM_VORTICITY_MOMENTUM_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_random_vorticity_momentum.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace random_vorticity_momentum_initial_condition {
    inline const std::string TYPE_NAME = "random_vorticity";
}

struct RandomVorticityMomentumInitialConditionCommand : public MomentumInitialConditionCommandBase {
    double n0 = 0.0;
    double sigma = 0.0;
    double urms = 0.0;
    int seed = 12345;

    void print(std::ostream& os) const override {
        const char* component_name = "?";
        if (direction == 0) {
            component_name = "x";
        } else if (direction == 1) {
            component_name = "y";
        }

        os << "  "
           << std::left << std::setw(25)
           << "Momentum" << ": Component " << component_name
           << " " << type
           << " n0 " << n0
           << " sigma " << sigma
           << " urms " << urms
           << " seed " << seed << '\n';
    }
};

class RandomVorticityMomentumInitialConditionStyle : public MomentumInitialConditionStyle {
private:
    const std::string name_ = random_vorticity_momentum_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MomentumInitialConditionCommandBase> parse_command(
        int direction,
        const InitialConditionArgs& args,
        const Params& /*params*/
    ) const override {
        auto cmd = std::make_shared<RandomVorticityMomentumInitialConditionCommand>();
        cmd->type = random_vorticity_momentum_initial_condition::TYPE_NAME;
        cmd->direction = direction;
        cmd->n0 = std::stod(args.get_required("n0"));
        cmd->sigma = std::stod(args.get_required("sigma"));
        cmd->urms = std::stod(args.get_required("urms"));
        cmd->seed = std::stoi(args.get_or("seed", "12345"));

        if (cmd->n0 <= 0.0) {
            throw std::runtime_error("set momentum random_vorticity requires positive n0.");
        }
        if (cmd->sigma <= 0.0) {
            throw std::runtime_error("set momentum random_vorticity requires positive sigma.");
        }
        if (cmd->urms <= 0.0) {
            throw std::runtime_error("set momentum random_vorticity requires positive urms.");
        }

        return cmd;
    }

    std::unique_ptr<MomentumInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    ) const override {
        const auto random_command = std::dynamic_pointer_cast<const RandomVorticityMomentumInitialConditionCommand>(command);
        if (!random_command) {
            throw std::runtime_error("RandomVorticityMomentumInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<RandomVorticityMomentumInitialCondition>(params, random_command);
    }
};

#endif
