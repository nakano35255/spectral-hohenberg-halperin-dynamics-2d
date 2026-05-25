#ifndef SFI_INITIAL_CONDITION_UNIFORM_ORDER_PARAMETER_STYLE_H
#define SFI_INITIAL_CONDITION_UNIFORM_ORDER_PARAMETER_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_uniform_order_parameter.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <string>

namespace uniform_order_parameter_initial_condition {
    inline const std::string TYPE_NAME = "uniform";
}

struct UniformOrderParameterInitialConditionCommand
    : public OrderParameterInitialConditionCommandBase {
    double value = 0.0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Order parameter" << ": Component " << component
           << " " << type
           << " value " << value << '\n';
    }
};

class UniformOrderParameterInitialConditionStyle : public OrderParameterInitialConditionStyle {
private:
    const std::string name_ = uniform_order_parameter_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<OrderParameterInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        (void)params;

        auto cmd = std::make_shared<UniformOrderParameterInitialConditionCommand>();
        cmd->value = std::stod(args.get_required("value"));
        cmd->type = uniform_order_parameter_initial_condition::TYPE_NAME;
        cmd->component = component;

        return cmd;
    }

    std::unique_ptr<OrderParameterInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    ) const override {
        const auto uniform_command = std::dynamic_pointer_cast<const UniformOrderParameterInitialConditionCommand>(command);
        if (!uniform_command) {
            throw std::runtime_error("UniformOrderParameterInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<UniformOrderParameterInitialCondition>(params, uniform_command);
    }
};

#endif
