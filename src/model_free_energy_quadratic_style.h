#ifndef SFI_MODEL_FREE_ENERGY_QUADRATIC_STYLE_H
#define SFI_MODEL_FREE_ENERGY_QUADRATIC_STYLE_H

#include "model_free_energy_quadratic.h"
#include "model_free_energy_registry.h"
#include "simulationinfo.h"

#include <cctype>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace quadratic_free_energy {
    inline const std::string TYPE_NAME = "quadratic";
}

struct QuadraticFreeEnergyCommand : FreeEnergyCommandBase {
    int num_order_parameters = 0;
    std::vector<double> coefficients;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Free Energy" << ": " << type;

        for (int q = 0; q < num_order_parameters; ++q) {
            os << " a[" << q << "] "
               << coefficients[static_cast<std::size_t>(q)];
        }

        os << '\n';
    }
};

class QuadraticFreeEnergyStyle : public FreeEnergyStyle {
private:
    const std::string name_ = quadratic_free_energy::TYPE_NAME;

    static bool parse_coefficient_key(const std::string& key, int& order_parameter) {
        if (key.size() < 4 || key.back() != ']') {
            return false;
        }
        if (key.rfind("a[", 0) != 0) {
            return false;
        }

        for (std::size_t i = 2; i < key.size() - 1; ++i) {
            if (!std::isdigit(static_cast<unsigned char>(key[i]))) {
                return false;
            }
        }

        order_parameter = std::stoi(key.substr(2, key.size() - 3));
        return true;
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<FreeEnergyCommandBase> create_default_command(const Params& params) const override {
        if (params.physics.num_order_parameters < 0) {
            throw std::runtime_error("model free_energy quadratic requires nonnegative number of order parameters.");
        }

        auto cmd = std::make_shared<QuadraticFreeEnergyCommand>();
        cmd->type = quadratic_free_energy::TYPE_NAME;
        cmd->num_order_parameters = params.physics.num_order_parameters;
        cmd->coefficients.assign(static_cast<std::size_t>(cmd->num_order_parameters), 0.0);
        return cmd;
    }

    void update_command(FreeEnergyCommandBase& command, const FreeEnergyArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<QuadraticFreeEnergyCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("QuadraticFreeEnergyStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            int order_parameter = -1;
            if (parse_coefficient_key(kv.first, order_parameter)) {
                if (order_parameter < 0 || order_parameter >= cmd->num_order_parameters) {
                    throw std::runtime_error("model free_energy quadratic coefficient index out of range: " + kv.first);
                }

                cmd->coefficients[static_cast<std::size_t>(order_parameter)] = std::stod(kv.second);
                continue;
            }

            throw std::runtime_error("unknown quadratic free energy coefficient: " + kv.first);
        }
    }

    std::unique_ptr<FreeEnergy> create(const Params& params, std::shared_ptr<const FreeEnergyCommandBase> command) const override {
        (void)params;

        const auto cmd = std::dynamic_pointer_cast<const QuadraticFreeEnergyCommand>(command);
        if (!cmd) {
            throw std::runtime_error("QuadraticFreeEnergyStyle received incompatible command.");
        }

        return std::make_unique<QuadraticFreeEnergy>(
            cmd->num_order_parameters,
            cmd->coefficients
        );
    }
};

#endif
