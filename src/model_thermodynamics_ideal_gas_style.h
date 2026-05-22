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

namespace ideal_gas_thermodynamics {
    inline const std::string TYPE_NAME = "ideal_gas";
}

struct IdealGasThermodynamicsCommand : ThermodynamicsCommandBase {
    double kB = 1.0;
    double temperature = 1.0;
    std::vector<double> mass;
    std::vector<double> reference_density;

    bool kB_set = false;
    bool temperature_set = false;
    std::vector<bool> mass_set;
    std::vector<bool> reference_density_set;

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

class IdealGasThermodynamicsStyle : public ThermodynamicsStyle {
private:
    const std::string name_ = ideal_gas_thermodynamics::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    void validate_command(const ThermodynamicsCommandBase& command, const ThermodynamicsRequirement& requirement, const Params& /*params*/) const override {
        const auto* cmd = dynamic_cast<const IdealGasThermodynamicsCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("IdealGasThermodynamicsStyle received incompatible command.");
        }

        const bool thermodynamics_required =
            requirement.free_energy ||
            requirement.chemical_potential ||
            requirement.bulk_pressure ||
            requirement.interfacial_force;

        if (!thermodynamics_required) {
            return;
        }

        if (requirement.interfacial_force) {
            throw std::runtime_error("model thermo ideal_gas does not support interfacial force.");
        }

        if (!cmd->kB_set) {
            throw std::runtime_error("model thermo ideal_gas requires kB for current fix settings.");
        }

        if (!cmd->temperature_set) {
            throw std::runtime_error("model thermo ideal_gas requires T for current fix settings.");
        }

        for (std::size_t k = 0; k < cmd->mass.size(); ++k) {
            if (!cmd->mass_set[k]) {
                throw std::runtime_error("model thermo ideal_gas requires m" + std::to_string(k) + " for current fix settings.");
            }
        }

        for (std::size_t k = 0; k < cmd->reference_density.size(); ++k) {
            if (!cmd->reference_density_set[k]) {
                throw std::runtime_error("model thermo ideal_gas requires rho" + std::to_string(k) + "_ref for current fix settings.");
            }
        }
    }

    std::shared_ptr<ThermodynamicsCommandBase> create_default_command(const Params& params) const override {
        if (params.physics.num_components <= 0) {
            throw std::runtime_error("model thermo ideal_gas requires components to be specified first.");
        }

        auto cmd = std::make_shared<IdealGasThermodynamicsCommand>();
        cmd->type = ideal_gas_thermodynamics::TYPE_NAME;

        const int n = params.physics.num_components;
        cmd->mass.assign(static_cast<std::size_t>(n), 1.0);
        cmd->reference_density.assign(static_cast<std::size_t>(n), 1.0);
        cmd->mass_set.assign(n, false);
        cmd->reference_density_set.assign(n, false);

        return cmd;
    }

    void update_command(ThermodynamicsCommandBase& command, const ThermodynamicsArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<IdealGasThermodynamicsCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("IdealGasThermodynamicsStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            const std::string& key = kv.first;
            const double value = std::stod(kv.second);

            if (key == "kB") {
                if (value <= 0.0) {
                    throw std::runtime_error("model thermo ideal_gas requires positive kB.");
                }
                cmd->kB = value;
                cmd->kB_set = true;
                continue;
            }

            if (key == "T") {
                if (value <= 0.0) {
                    throw std::runtime_error("model thermo ideal_gas requires positive T.");
                }
                cmd->temperature = value;
                cmd->temperature_set = true;
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
                    cmd->mass_set[k] = true;
                    matched = true;
                    break;
                }

                const std::string rho_ref_key = "rho" + std::to_string(k) + "_ref";
                if (key == rho_ref_key) {
                    if (value <= 0.0) {
                        throw std::runtime_error("model thermo ideal_gas requires positive " + key + ".");
                    }
                    cmd->reference_density[k] = value;
                    cmd->reference_density_set[k] = true;
                    matched = true;
                    break;
                }
            }

            if (!matched) {
                throw std::runtime_error("unknown ideal_gas thermodynamics coefficient: " + key);
            }
        }
    }

    std::unique_ptr<Thermodynamics> create(const Params& params, std::shared_ptr<const ThermodynamicsCommandBase> command) const override {
        (void)params;

        const auto cmd = std::dynamic_pointer_cast<const IdealGasThermodynamicsCommand>(command);
        if (!cmd) {
            throw std::runtime_error("IdealGasThermodynamicsStyle received incompatible command.");
        }

        return std::make_unique<IdealGasThermodynamics>(
            cmd->kB,
            cmd->temperature,
            cmd->mass,
            cmd->reference_density
        );
    }
};

#endif
