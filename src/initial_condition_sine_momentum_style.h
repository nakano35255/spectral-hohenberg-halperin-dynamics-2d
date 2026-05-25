#ifndef SFI_INITIAL_CONDITION_SINE_MOMENTUM_STYLE_H
#define SFI_INITIAL_CONDITION_SINE_MOMENTUM_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_sine_momentum.h"

#include <cmath>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace sine_momentum_initial_condition {
    inline const std::string TYPE_NAME = "sine";
}

struct SineMomentumInitialConditionCommand : public MomentumInitialConditionCommandBase {
    double base = 0.0;
    double amplitude = 0.0;
    int nkx = 0;
    int nky = 0;

    void print(std::ostream& os) const override {
        const char* direction_name = "?";
        if (direction == 0) {
            direction_name = "x";
        } else if (direction == 1) {
            direction_name = "y";
        }

        os << "  "
           << std::left << std::setw(25)
           << "Momentum" << ": Direction " << direction_name
           << " " << type
           << " base " << base
           << " amplitude " << amplitude
           << " nkx " << nkx
           << " nky " << nky << '\n';
    }
};

class SineMomentumInitialConditionStyle : public MomentumInitialConditionStyle {
private:
    const std::string name_ = sine_momentum_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MomentumInitialConditionCommandBase> parse_command(
        int direction,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<SineMomentumInitialConditionCommand>();
        cmd->base = std::stod(args.get_required("base"));
        cmd->amplitude = std::stod(args.get_required("amplitude"));
        cmd->nkx = std::stoi(args.get_required("nkx"));
        cmd->nky = std::stoi(args.get_or("nky", "0"));
        cmd->type = sine_momentum_initial_condition::TYPE_NAME;
        cmd->direction = direction;

        if (cmd->nkx < 0 || cmd->nkx > params.grid.active_num_grid[0] / 2) {
            throw std::runtime_error("set momentum sine requires 0 <= nkx <= Nx/2");
        }
        if (cmd->nky < -params.grid.active_num_grid[1] / 2 || cmd->nky > params.grid.active_num_grid[1] / 2) {
            throw std::runtime_error("set momentum sine requires -Ny/2 <= nky <= Ny/2");
        }
        if (cmd->nkx == 0 && cmd->nky == 0) {
            throw std::runtime_error("set momentum sine requires a nonzero wave number");
        }

        const bool self_conjugate_kx = (cmd->nkx == 0 || cmd->nkx == params.grid.active_num_grid[0] / 2);
        const bool self_conjugate_ky = (cmd->nky == 0 || std::abs(cmd->nky) == params.grid.active_num_grid[1] / 2);
        if (self_conjugate_kx && self_conjugate_ky) {
            throw std::runtime_error("set momentum sine mode is zero on the grid; choose a non-self-conjugate wave number");
        }

        return cmd;
    }

    std::unique_ptr<MomentumInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    ) const override {
        const auto sine_command = std::dynamic_pointer_cast<const SineMomentumInitialConditionCommand>(command);
        if (!sine_command) {
            throw std::runtime_error("SineMomentumInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<SineMomentumInitialCondition>(params, sine_command);
    }
};

#endif
