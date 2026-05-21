#ifndef SFI_MODEL_THERMODYNAMICS_IDEAL_GAS_STYLE_H
#define SFI_MODEL_THERMODYNAMICS_IDEAL_GAS_STYLE_H

#include "model_thermodynamics_ideal_gas.h"
#include "model_thermodynamics_registry.h"
#include "simulationinfo.h"

#include <cstddef>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ideal_gas_thermodynamics_model {
    inline const std::string TYPE_NAME = "ideal_gas";
}

struct IdealGasThermodynamicsModelCommand : ThermodynamicsModelCommandBase {
    double kB = 1.0;
    double temperature = 1.0;
    std::vector<double> mass;
    std::vector<double> reference_density;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Thermodynamics" << ": " << type
           << " kB " << kB
           << " T " << temperature;

        for (std::size_t k = 0; k < mass.size(); ++k) {
            os << " m" << k << " " << mass[k]
               << " rho" << k << "_ref " << reference_density[k];
        }

        os << '\n';
    }
};

class IdealGasThermodynamicsModelStyle : public ThermodynamicsModelStyle {
private:
    const std::string name_ = ideal_gas_thermodynamics_model::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<ThermodynamicsModelCommandBase> create_default_command(
        const Params& params
    ) const override {
        if (params.physics.num_components <= 0) {
            throw std::runtime_error("model thermo ideal_gas requires components to be specified first.");
        }

        auto cmd = std::make_shared<IdealGasThermodynamicsModelCommand>();
        cmd->type = ideal_gas_thermodynamics_model::TYPE_NAME;

        const int n = params.physics.num_components;
        cmd->mass.assign(static_cast<std::size_t>(n), 1.0);
        cmd->reference_density.assign(static_cast<std::size_t>(n), 1.0);

        return cmd;
    }

    void update_command(
        ThermodynamicsModelCommandBase& command,
        const ThermodynamicsModelArgs& args,
        const Params& /*params*/
    ) const override {
        auto* cmd = dynamic_cast<IdealGasThermodynamicsModelCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("IdealGasThermodynamicsModelStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            const std::string& key = kv.first;
            const double value = std::stod(kv.second);

            if (key == "kB") {
                if (value <= 0.0) {
                    throw std::runtime_error("model thermo ideal_gas requires positive kB.");
                }
                cmd->kB = value;
                continue;
            }

            if (key == "T") {
                if (value <= 0.0) {
                    throw std::runtime_error("model thermo ideal_gas requires positive T.");
                }
                cmd->temperature = value;
                continue;
            }

            bool matched = false;
            for (std::size_t k = 0; k < cmd->mass.size(); ++k) {
                const std::string mass_key = "m" + std::to_string(k);
                if (key == mass_key) {
                    if (value <= 0.0) {
                        throw std::runtime_error("model thermo ideal_gas requires positive " + key + ".");
                    }
                    cmd->mass[k] = value;
                    matched = true;
                    break;
                }

                const std::string rho_ref_key = "rho" + std::to_string(k) + "_ref";
                if (key == rho_ref_key) {
                    if (value <= 0.0) {
                        throw std::runtime_error("model thermo ideal_gas requires positive " + key + ".");
                    }
                    cmd->reference_density[k] = value;
                    matched = true;
                    break;
                }
            }

            if (!matched) {
                throw std::runtime_error("unknown ideal_gas thermodynamics coefficient: " + key);
            }
        }
    }

    std::unique_ptr<ThermodynamicsModel> create_model(
        const Params& params,
        std::shared_ptr<const ThermodynamicsModelCommandBase> command
    ) const override {
        (void)params;

        const auto cmd =
            std::dynamic_pointer_cast<const IdealGasThermodynamicsModelCommand>(command);
        if (!cmd) {
            throw std::runtime_error("IdealGasThermodynamicsModelStyle received incompatible command.");
        }

        return std::make_unique<IdealGasThermodynamicsModel>(
            cmd->kB,
            cmd->temperature,
            cmd->mass,
            cmd->reference_density
        );
    }
};

#endif
