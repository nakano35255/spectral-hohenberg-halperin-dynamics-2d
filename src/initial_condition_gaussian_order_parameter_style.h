#ifndef SFI_INITIAL_CONDITION_GAUSSIAN_ORDER_PARAMETER_STYLE_H
#define SFI_INITIAL_CONDITION_GAUSSIAN_ORDER_PARAMETER_STYLE_H

#include "initial_condition_gaussian_order_parameter.h"
#include "initial_condition_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace gaussian_order_parameter_initial_condition {
    inline const std::string TYPE_NAME = "gaussian";
}

struct GaussianOrderParameterInitialConditionCommand : public OrderParameterInitialConditionCommandBase {
    double base = 0.0;
    double amplitude = 0.0;
    double x0 = 0.0;
    double y0 = 0.0;
    double sigma_x = 1.0;
    double sigma_y = 1.0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Order parameter" << ": Component " << component
           << " " << type
           << " base " << base
           << " amplitude " << amplitude
           << " x0 " << x0
           << " y0 " << y0
           << " sigma_x " << sigma_x
           << " sigma_y " << sigma_y << '\n';
    }
};

class GaussianOrderParameterInitialConditionStyle : public OrderParameterInitialConditionStyle {
private:
    const std::string name_ = gaussian_order_parameter_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<OrderParameterInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<GaussianOrderParameterInitialConditionCommand>();
        cmd->base = std::stod(args.get_required("base"));
        cmd->amplitude = std::stod(args.get_required("amplitude"));
        cmd->x0 = std::stod(args.get_required("x0"));
        cmd->y0 = std::stod(args.get_required("y0"));

        if (args.has("sigma")) {
            const double sigma = std::stod(args.get_required("sigma"));
            cmd->sigma_x = sigma;
            cmd->sigma_y = sigma;
        } else {
            cmd->sigma_x = std::stod(args.get_required("sigma_x"));
            cmd->sigma_y = std::stod(args.get_required("sigma_y"));
        }

        cmd->type = gaussian_order_parameter_initial_condition::TYPE_NAME;
        cmd->component = component;

        if (cmd->sigma_x <= 0.0 || cmd->sigma_y <= 0.0) {
            throw std::runtime_error("set order_parameter gaussian requires positive sigma values");
        }
        if (params.grid.length[0] <= 0.0 || params.grid.length[1] <= 0.0) {
            throw std::runtime_error("set order_parameter gaussian requires positive domain lengths");
        }

        return cmd;
    }

    std::unique_ptr<OrderParameterInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
    ) const override {
        const auto gaussian_command = std::dynamic_pointer_cast<const GaussianOrderParameterInitialConditionCommand>(command);
        if (!gaussian_command) {
            throw std::runtime_error("GaussianOrderParameterInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<GaussianOrderParameterInitialCondition>(params, gaussian_command);
    }
};

#endif
