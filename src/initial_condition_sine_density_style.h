#ifndef SFI_INITIAL_CONDITION_SINE_DENSITY_STYLE_H
#define SFI_INITIAL_CONDITION_SINE_DENSITY_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_sine_density.h"

#include <cmath>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace sine_density_initial_condition {
    inline const std::string TYPE_NAME = "sine";
}

struct SineDensityInitialConditionCommand : public DensityInitialConditionCommandBase {
    double base = 0.0;
    double amplitude = 0.0;
    int nkx = 0;
    int nky = 0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Density"
           << " " << type
           << " base " << base
           << " amplitude " << amplitude
           << " nkx " << nkx
           << " nky " << nky << '\n';
    }
};

class SineDensityInitialConditionStyle : public DensityInitialConditionStyle {
private:
    const std::string name_ = sine_density_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<DensityInitialConditionCommandBase> parse_command(
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<SineDensityInitialConditionCommand>();
        cmd->base = std::stod(args.get_required("base"));
        cmd->amplitude = std::stod(args.get_required("amplitude"));
        cmd->nkx = std::stoi(args.get_required("nkx"));
        cmd->nky = std::stoi(args.get_or("nky", "0"));
        cmd->type = sine_density_initial_condition::TYPE_NAME;

        if (cmd->base - std::abs(cmd->amplitude) < 0.0) {
            throw std::runtime_error("set density sine requires base >= abs(amplitude)");
        }
        if (cmd->nkx < 0 || cmd->nkx >= params.grid.active_num_grid[0] / 2) {
            throw std::runtime_error("set density sine requires 0 <= nkx < Nx/2; Nyquist line is fixed to zero");
        }
        if (cmd->nky <= -params.grid.active_num_grid[1] / 2 || cmd->nky >=  params.grid.active_num_grid[1] / 2) {
            throw std::runtime_error("set density sine requires -Ny/2 < nky < Ny/2; Nyquist line is fixed to zero");
        }
        if (cmd->nkx == 0 && cmd->nky == 0) {
            throw std::runtime_error("set density sine requires a nonzero wave number");
        }

        return cmd;
    }

    std::unique_ptr<DensityInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    ) const override {
        const auto sine_command = std::dynamic_pointer_cast<const SineDensityInitialConditionCommand>(command);
        if (!sine_command) {
            throw std::runtime_error("SineDensityInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<SineDensityInitialCondition>(params, sine_command);
    }
};

#endif
