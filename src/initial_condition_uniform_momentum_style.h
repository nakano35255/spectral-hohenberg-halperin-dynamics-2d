#ifndef SHHD_INITIAL_CONDITION_UNIFORM_MOMENTUM_STYLE_H
#define SHHD_INITIAL_CONDITION_UNIFORM_MOMENTUM_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_uniform_momentum.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace uniform_momentum_initial_condition {
    inline const std::string TYPE_NAME = "uniform";
}

struct UniformMomentumInitialConditionCommand : public MomentumInitialConditionCommandBase {
    double value = 0.0;

    UniformMomentumInitialConditionCommand(int direction_in, double value_in)
        : value(value_in) {
        type = uniform_momentum_initial_condition::TYPE_NAME;
        direction = direction_in;
    }

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
           << " value " << value << '\n';
    }
};

class UniformMomentumInitialConditionStyle : public MomentumInitialConditionStyle {
private:
    const std::string name_ = uniform_momentum_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override{
        return name_;
    };

    std::shared_ptr<MomentumInitialConditionCommandBase> parse_command(int component, const InitialConditionArgs& args, const Params& params) const override {
        (void)params;

        const std::string value_text = args.get_required("value");
        const double value = std::stod(value_text);

        return std::make_shared<UniformMomentumInitialConditionCommand>(component, value);
    };

    std::unique_ptr<MomentumInitialCondition> create_initial_condition(const Params& params, std::shared_ptr<const MomentumInitialConditionCommandBase> command) const override {
        (void)params;

        const auto uniform_command = std::dynamic_pointer_cast<const UniformMomentumInitialConditionCommand>(command);

        if (!uniform_command) {
            throw std::runtime_error("UniformMomentumInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<UniformMomentumInitialCondition>(params, uniform_command);
    }
};

#endif
