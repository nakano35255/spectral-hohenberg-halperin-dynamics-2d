#ifndef SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_STYLE_H
#define SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_STYLE_H

#include "model_transport_coefficient_constant.h"
#include "model_transport_coefficient_registry.h"
#include "simulationinfo.h"

#include <cctype>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace constant_transport_coefficient {
    inline const std::string TYPE_NAME = "constant";
}

struct ConstantTransportCoefficientCommand : TransportCoefficientCommandBase {
    int num_components = 0;
    double shear_viscosity = 0.0;
    double bulk_viscosity = 0.0;
    std::vector<double> mobility;

    bool shear_viscosity_set = false;
    bool bulk_viscosity_set = false;
    std::vector<bool> mobility_set;

    void print(std::ostream& os) const override {
        os << "  "
           << std::left << std::setw(25)
           << "Transport Coefficients" << ": " << type
           << " eta " << shear_viscosity
           << " zeta " << bulk_viscosity;

        for (int i = 0; i < num_components; ++i) {
            for (int j = 0; j < num_components; ++j) {
                const std::size_t index = static_cast<std::size_t>(i) * static_cast<std::size_t>(num_components) + static_cast<std::size_t>(j);
                os << " L" << i << j << " " << mobility[index];
            }
        }

        os << '\n';
    }
};

class ConstantTransportCoefficientStyle : public TransportCoefficientStyle {
private:
    const std::string name_ = constant_transport_coefficient::TYPE_NAME;

    static std::size_t mobility_index(int row, int col, int num_components) {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(num_components) + static_cast<std::size_t>(col);
    }

    static bool parse_mobility_key(const std::string& key, int& row, int& col) {
        if (key.size() != 3 || key[0] != 'L') {
            return false;
        }

        if (!std::isdigit(static_cast<unsigned char>(key[1])) || !std::isdigit(static_cast<unsigned char>(key[2]))) {
            return false;
        }

        row = key[1] - '0';
        col = key[2] - '0';
        return true;
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<TransportCoefficientCommandBase> create_default_command(const Params& params) const override {
        if (params.physics.num_components <= 0) {
            throw std::runtime_error("model transport constant requires components to be specified first.");
        }

        auto cmd = std::make_shared<ConstantTransportCoefficientCommand>();
        cmd->type = constant_transport_coefficient::TYPE_NAME;
        cmd->num_components = params.physics.num_components;

        const std::size_t matrix_size = static_cast<std::size_t>(cmd->num_components) * static_cast<std::size_t>(cmd->num_components);
        cmd->mobility.assign(matrix_size, 0.0);
        cmd->mobility_set.assign(matrix_size, false);

        return cmd;
    }

    void update_command(TransportCoefficientCommandBase& command, const TransportCoefficientArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<ConstantTransportCoefficientCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("ConstantTransportCoefficientStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            const std::string& key = kv.first;
            const double value = std::stod(kv.second);

            if (key == "eta" || key == "shear_viscosity") {
                if (value < 0.0) {
                    throw std::runtime_error("model transport constant requires nonnegative eta.");
                }
                cmd->shear_viscosity = value;
                cmd->shear_viscosity_set = true;
                continue;
            }

            if (key == "zeta" || key == "bulk_viscosity") {
                if (value < 0.0) {
                    throw std::runtime_error("model transport constant requires nonnegative zeta.");
                }
                cmd->bulk_viscosity = value;
                cmd->bulk_viscosity_set = true;
                continue;
            }

            int row = -1;
            int col = -1;
            if (parse_mobility_key(key, row, col)) {
                if (row < 0 || row >= cmd->num_components || col < 0 || col >= cmd->num_components) {
                    throw std::runtime_error("model transport constant mobility index out of range: " + key);
                }

                cmd->mobility[mobility_index(row, col, cmd->num_components)] = value;
                cmd->mobility_set[mobility_index(row, col, cmd->num_components)] = true;
                continue;
            }

            throw std::runtime_error("unknown constant transport coefficient: " + key);
        }
    }

    std::unique_ptr<TransportCoefficient> create(const Params& params, std::shared_ptr<const TransportCoefficientCommandBase> command) const override {
        (void)params;

        const auto cmd = std::dynamic_pointer_cast<const ConstantTransportCoefficientCommand>(command);
        if (!cmd) {
            throw std::runtime_error("ConstantTransportCoefficientStyle received incompatible command.");
        }

        return std::make_unique<ConstantTransportCoefficient>(
            cmd->num_components,
            cmd->shear_viscosity,
            cmd->bulk_viscosity,
            cmd->mobility
        );
    }
};

#endif
