#ifndef SHHD_MODEL_THERMODYNAMICS_LINEAR_EOS_STYLE_H
#define SHHD_MODEL_THERMODYNAMICS_LINEAR_EOS_STYLE_H

#include "model_thermodynamics_linear_eos.h"
#include "model_thermodynamics_registry.h"
#include "simulationinfo.h"

#include <memory>
#include <string>

namespace linear_eos_thermodynamics {
    inline const std::string TYPE_NAME = "linear_eos";
}

struct LinearEOSThermodynamicsCommand : ThermodynamicsCommandBase {
    void print(std::ostream& os) const override {
        print_common(os);
    }
};

class LinearEOSThermodynamicsStyle : public ThermodynamicsStyle {
private:
    const std::string name_ = linear_eos_thermodynamics::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<ThermodynamicsCommandBase> create_default_command(const Params& /*params*/) const override {
        auto cmd = std::make_shared<LinearEOSThermodynamicsCommand>();
        cmd->type = linear_eos_thermodynamics::TYPE_NAME;
        return cmd;
    }

    void update_command(ThermodynamicsCommandBase& command, const ThermodynamicsArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<LinearEOSThermodynamicsCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("LinearEOSThermodynamicsStyle received incompatible command.");
        }
        for (const auto& kv : args.entries) {
            const std::string& key = kv.first;

            if (update_common_command(*cmd, key, kv.second)) {
                continue;
            }

            throw std::runtime_error("unknown linear_eos thermodynamics coefficient: " + key);
        }
    }

    std::unique_ptr<Thermodynamics> create(const Params& params, std::shared_ptr<const ThermodynamicsCommandBase> command) const override {
        (void)params;

        const auto cmd = std::dynamic_pointer_cast<const LinearEOSThermodynamicsCommand>(command);
        if (!cmd) {
            throw std::runtime_error("LinearEOSThermodynamicsStyle received incompatible command.");
        }

        return std::make_unique<LinearEOSThermodynamics>(cmd->sound_speed);
    }
};

#endif
