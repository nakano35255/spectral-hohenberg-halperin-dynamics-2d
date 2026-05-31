#ifndef SHHD_INITIAL_CONDITION_SINE_ORDER_PARAMETER_STYLE_H
#define SHHD_INITIAL_CONDITION_SINE_ORDER_PARAMETER_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_sine_order_parameter.h"

#include <cmath>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace sine_order_parameter_initial_condition {
    inline const std::string TYPE_NAME = "sine";
}

struct SineOrderParameterInitialConditionCommand : public OrderParameterInitialConditionCommandBase {
    double base = 0.0;
    double amplitude = 0.0;
    int nkx = 0;
    int nky = 0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Order parameter" << ": Component " << component
           << " " << type
           << " base " << base
           << " amplitude " << amplitude
           << " nkx " << nkx
           << " nky " << nky << '\n';
    }
};

class SineOrderParameterInitialConditionStyle : public OrderParameterInitialConditionStyle {
private:
    const std::string name_ = sine_order_parameter_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<OrderParameterInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<SineOrderParameterInitialConditionCommand>();
        cmd->base = std::stod(args.get_or("base", "0.0"));
        cmd->amplitude = std::stod(args.get_required("amplitude"));
        cmd->nkx = std::stoi(args.get_required("nkx"));
        cmd->nky = std::stoi(args.get_or("nky", "0"));
        cmd->type = sine_order_parameter_initial_condition::TYPE_NAME;
        cmd->component = component;

        if (cmd->nkx < 0 || cmd->nkx >= params.grid.active_num_grid[0] / 2) {
            throw std::runtime_error("set order_parameter sine requires 0 <= nkx < Nx/2; Nyquist line is fixed to zero");
        }

        if (cmd->nky <= -params.grid.active_num_grid[1] / 2 || cmd->nky >=  params.grid.active_num_grid[1] / 2) {
            throw std::runtime_error("set order_parameter sine requires -Ny/2 < nky < Ny/2; Nyquist line is fixed to zero");
        }
        if (cmd->nkx == 0 && cmd->nky == 0) {
            throw std::runtime_error("set order_parameter sine requires a nonzero wave number");
        }
        
        return cmd;
    }

    std::unique_ptr<OrderParameterInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    ) const override {
        const auto sine_command = std::dynamic_pointer_cast<const SineOrderParameterInitialConditionCommand>(command);
        if (!sine_command) {
            throw std::runtime_error("SineOrderParameterInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<SineOrderParameterInitialCondition>(params, sine_command);
    }
};

#endif
