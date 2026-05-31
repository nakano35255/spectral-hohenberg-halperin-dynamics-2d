#ifndef SHHD_INITIAL_CONDITION_UNIFORM_DENSITY_STYLE_H
#define SHHD_INITIAL_CONDITION_UNIFORM_DENSITY_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_uniform_density.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <string>

namespace uniform_density_initial_condition {
    inline const std::string TYPE_NAME = "uniform";
}

struct UniformDensityInitialConditionCommand : public DensityInitialConditionCommandBase {
    double value = 0.0;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Density"
           << " " << type
           << " value " << value << '\n';
    }
};

class UniformDensityInitialConditionStyle : public DensityInitialConditionStyle {
private:
    const std::string name_ = uniform_density_initial_condition::TYPE_NAME;

public:
    const std::string& type_name() const override{
        return name_;
    };

    std::shared_ptr<DensityInitialConditionCommandBase> parse_command(const InitialConditionArgs& args, const Params& /*params*/) const override {
        auto cmd = std::make_shared<UniformDensityInitialConditionCommand>();
        cmd->value = std::stod(args.get_required("value"));
        cmd->type = uniform_density_initial_condition::TYPE_NAME;
        
        if (cmd->value < 0.0) {
            throw std::runtime_error("set density uniform requires nonnegative value");
        }

        return cmd;
    };

    std::unique_ptr<DensityInitialCondition> create_initial_condition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command) const override {
        const auto uniform_command = std::dynamic_pointer_cast<const UniformDensityInitialConditionCommand>(command);
        if (!uniform_command) {
            throw std::runtime_error("UniformDensityInitialConditionStyle received an incompatible command.");
        }

        return std::make_unique<UniformDensityInitialCondition>(params, uniform_command);
    }
};

#endif
