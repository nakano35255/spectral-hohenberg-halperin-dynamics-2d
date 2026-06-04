#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_ORDER_PARAMETER_STYLE_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_ORDER_PARAMETER_STYLE_H

#include "initial_condition_equilibrium_gaussian_order_parameter.h"
#include "initial_condition_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace equilibrium_gaussian_order_parameter_initial_condition {
    inline const std::string TYPE_NAME = "equilibrium/gaussian";
}

struct EquilibriumGaussianOrderParameterInitialConditionCommand : public OrderParameterInitialConditionCommandBase {
    int seed = 12345;
    double kBT = 1.0;
    double psi0 = 0.0;
    double k0_coefficient = 0.0;
    double k2_coefficient = 0.0;
    double k4_coefficient = 0.0;

    void print(std::ostream& os) const override {
        os << "  " << std::left << std::setw(25)
           << "Order parameter" << ": Component " << component
           << " " << type
           << " psi0 " << psi0
           << " kBT " << kBT
           << " k0 " << k0_coefficient
           << " k2 " << k2_coefficient
           << " k4 " << k4_coefficient
           << " seed " << seed << '\n';
    }
};

class EquilibriumGaussianOrderParameterInitialConditionStyle : public OrderParameterInitialConditionStyle {
private:
    const std::string name_ = equilibrium_gaussian_order_parameter_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override { return name_; }

    std::shared_ptr<OrderParameterInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        (void)params;

        auto cmd = std::make_shared<EquilibriumGaussianOrderParameterInitialConditionCommand>();
        cmd->type = name_;
        cmd->component = component;
        cmd->seed = std::stoi(args.get_or("seed", "12345"));
        cmd->kBT = std::stod(args.get_required("kBT"));
        cmd->psi0 = std::stod(args.get_or("psi0", "0.0"));
        cmd->k0_coefficient = std::stod(args.get_required("k0"));
        cmd->k2_coefficient = std::stod(args.get_or("k2", "0.0"));
        cmd->k4_coefficient = std::stod(args.get_or("k4", "0.0"));

        if (cmd->kBT < 0.0) {
            throw std::runtime_error("set order_parameter " + name_ + " requires kBT >= 0.");
        }

        return cmd;
    }

    std::unique_ptr<OrderParameterInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    ) const override {
        return std::make_unique<EquilibriumGaussianOrderParameterInitialCondition>(params, command);
    }
};

#endif
