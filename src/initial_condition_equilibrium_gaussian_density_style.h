#ifndef SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_DENSITY_STYLE_H
#define SHHD_INITIAL_CONDITION_EQUILIBRIUM_GAUSSIAN_DENSITY_STYLE_H

#include "initial_condition_equilibrium_gaussian_density.h"
#include "initial_condition_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace equilibrium_gaussian_density_initial_condition {
    inline const std::string TYPE_NAME = "equilibrium/gaussian";
}

struct EquilibriumGaussianDensityInitialConditionCommand : public DensityInitialConditionCommandBase {
    int seed = 12345;
    double kBT = 1.0;
    double rho0 = 1.0;
    double cT = 1.0;
    double pressure_coefficient = 0.0;

    void print(std::ostream& os) const override {
        os << "  " << std::left << std::setw(25)
           << "Density" << ": " << type
           << " rho0 " << rho0
           << " kBT " << kBT
           << " cT " << cT
           << " seed " << seed << '\n';
    }
};

class EquilibriumGaussianDensityInitialConditionStyle : public DensityInitialConditionStyle {
private:
    const std::string name_ = equilibrium_gaussian_density_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override { return name_; }

    std::shared_ptr<DensityInitialConditionCommandBase> parse_command(
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        (void)params;

        auto cmd = std::make_shared<EquilibriumGaussianDensityInitialConditionCommand>();
        cmd->type = name_;
        cmd->seed = std::stoi(args.get_or("seed", "12345"));
        cmd->kBT = std::stod(args.get_required("kBT"));
        cmd->rho0 = std::stod(args.get_required("rho0"));
        cmd->cT = std::stod(args.get_required("cT"));
        cmd->pressure_coefficient = cmd->cT * cmd->cT;

        if (cmd->kBT < 0.0) {
            throw std::runtime_error("set density " + name_ + " requires kBT >= 0.");
        }
        if (cmd->rho0 <= 0.0) {
            throw std::runtime_error("set density " + name_ + " requires rho0 > 0.");
        }
        if (cmd->cT <= 0.0) {
            throw std::runtime_error("set density " + name_ + " requires cT > 0.");
        }

        return cmd;
    }

    std::unique_ptr<DensityInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    ) const override {
        return std::make_unique<EquilibriumGaussianDensityInitialCondition>(params, command);
    }
};

#endif
