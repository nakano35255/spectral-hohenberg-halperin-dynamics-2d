#ifndef SHHD_MODEL_FREE_ENERGY_PHI4_STYLE_H
#define SHHD_MODEL_FREE_ENERGY_PHI4_STYLE_H

#include "model_free_energy_phi4.h"
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

namespace phi4_free_energy {
    inline const std::string TYPE_NAME = "phi4";
}

struct Phi4FreeEnergyCommand : FreeEnergyCommandBase {
    int num_order_parameters = 0;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> u;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Free Energy" << ": " << type;

        for (int q = 0; q < num_order_parameters; ++q) {
            const std::size_t index = static_cast<std::size_t>(q);
            os << " a[" << q << "] " << a[index]
               << " b[" << q << "] " << b[index]
               << " u[" << q << "] " << u[index];
        }

        os << '\n';
    }
};

class Phi4FreeEnergyStyle : public FreeEnergyStyle {
private:
    const std::string name_ = phi4_free_energy::TYPE_NAME;

    static bool parse_indexed_key(const std::string& key, const std::string& prefix, int& order_parameter) {
        const std::string head = prefix + "[";
        if (key.size() <= head.size() + 1 || key.back() != ']') {
            return false;
        }
        if (key.rfind(head, 0) != 0) {
            return false;
        }

        for (std::size_t i = head.size(); i < key.size() - 1; ++i) {
            if (!std::isdigit(static_cast<unsigned char>(key[i]))) {
                return false;
            }
        }

        order_parameter = std::stoi(key.substr(head.size(), key.size() - head.size() - 1));
        return true;
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<FreeEnergyCommandBase> create_default_command(const Params& params) const override {
        if (params.physics.num_order_parameters < 0) {
            throw std::runtime_error("model free_energy phi4 requires nonnegative number of order parameters.");
        }

        auto cmd = std::make_shared<Phi4FreeEnergyCommand>();
        cmd->type = phi4_free_energy::TYPE_NAME;
        cmd->num_order_parameters = params.physics.num_order_parameters;
        cmd->a.assign(static_cast<std::size_t>(cmd->num_order_parameters), 0.0);
        cmd->b.assign(static_cast<std::size_t>(cmd->num_order_parameters), 0.0);
        cmd->u.assign(static_cast<std::size_t>(cmd->num_order_parameters), 0.0);
        return cmd;
    }

    void update_command(FreeEnergyCommandBase& command, const FreeEnergyArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<Phi4FreeEnergyCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("Phi4FreeEnergyStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            int order_parameter = -1;
            const double value = std::stod(kv.second);

            if (parse_indexed_key(kv.first, "a", order_parameter)) {
                if (order_parameter < 0 || order_parameter >= cmd->num_order_parameters) {
                    throw std::runtime_error("model free_energy phi4 coefficient index out of range: " + kv.first);
                }
                cmd->a[static_cast<std::size_t>(order_parameter)] = value;
                continue;
            }

            if (parse_indexed_key(kv.first, "b", order_parameter)) {
                if (order_parameter < 0 || order_parameter >= cmd->num_order_parameters) {
                    throw std::runtime_error("model free_energy phi4 coefficient index out of range: " + kv.first);
                }
                cmd->b[static_cast<std::size_t>(order_parameter)] = value;
                continue;
            }

            if (parse_indexed_key(kv.first, "u", order_parameter)) {
                if (order_parameter < 0 || order_parameter >= cmd->num_order_parameters) {
                    throw std::runtime_error("model free_energy phi4 coefficient index out of range: " + kv.first);
                }
                if (value < 0.0) {
                    throw std::runtime_error("model free_energy phi4 requires nonnegative " + kv.first + ".");
                }
                cmd->u[static_cast<std::size_t>(order_parameter)] = value;
                continue;
            }

            throw std::runtime_error("unknown phi4 free energy coefficient: " + kv.first);
        }
    }

    std::unique_ptr<FreeEnergy> create(const Params& params, std::shared_ptr<const FreeEnergyCommandBase> command) const override {
        (void)params;

        const auto cmd = std::dynamic_pointer_cast<const Phi4FreeEnergyCommand>(command);
        if (!cmd) {
            throw std::runtime_error("Phi4FreeEnergyStyle received incompatible command.");
        }

        return std::make_unique<Phi4FreeEnergy>(
            cmd->num_order_parameters,
            cmd->a,
            cmd->b,
            cmd->u
        );
    }
};

#endif
