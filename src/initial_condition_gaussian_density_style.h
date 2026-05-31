#ifndef SHHD_INITIAL_CONDITION_GAUSSIAN_DENSITY_STYLE_H
#define SHHD_INITIAL_CONDITION_GAUSSIAN_DENSITY_STYLE_H

#include "initial_condition_gaussian_density.h"
#include "initial_condition_registry.h"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace gaussian_density_initial_condition {
    inline const std::string TYPE_NAME = "gaussian";
}

struct GaussianDensityInitialConditionCommand : public DensityInitialConditionCommandBase {
    double base = 0.0;
    double amplitude = 0.0;
    double x0 = 0.0;
    double y0 = 0.0;
    double sigma_x = 1.0;
    double sigma_y = 1.0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Density"
           << " " << type
           << " base " << base
           << " amplitude " << amplitude
           << " x0 " << x0
           << " y0 " << y0
           << " sigma_x " << sigma_x
           << " sigma_y " << sigma_y << '\n';
    }
};

class GaussianDensityInitialConditionStyle : public DensityInitialConditionStyle {
private:
    const std::string name_ = gaussian_density_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<DensityInitialConditionCommandBase> parse_command(
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<GaussianDensityInitialConditionCommand>();
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

        cmd->type = gaussian_density_initial_condition::TYPE_NAME;

        if (cmd->sigma_x <= 0.0 || cmd->sigma_y <= 0.0) {
            throw std::runtime_error("set density gaussian requires positive sigma values");
        }
        if (cmd->base + std::min(0.0, cmd->amplitude) < 0.0) {
            throw std::runtime_error("set density gaussian requires nonnegative density everywhere");
        }
        if (params.grid.length[0] <= 0.0 || params.grid.length[1] <= 0.0) {
            throw std::runtime_error("set density gaussian requires positive domain lengths");
        }

        return cmd;
    }

    std::unique_ptr<DensityInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    ) const override {
        const auto gaussian_command = std::dynamic_pointer_cast<const GaussianDensityInitialConditionCommand>(command);
        if (!gaussian_command) {
            throw std::runtime_error("GaussianDensityInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<GaussianDensityInitialCondition>(params, gaussian_command);
    }
};

#endif
