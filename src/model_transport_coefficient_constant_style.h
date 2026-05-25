#ifndef SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_STYLE_H
#define SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_STYLE_H

#include "model_transport_coefficient_constant.h"
#include "model_transport_coefficient_registry.h"
#include "simulationinfo.h"

#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace constant_transport_coefficient {
    inline const std::string TYPE_NAME = "constant";
}

struct ConstantTransportCoefficientCommand : TransportCoefficientCommandBase {
    void print(std::ostream& os) const override {
        print_common(os);
    }
};

class ConstantTransportCoefficientStyle : public TransportCoefficientStyle {
private:
    const std::string name_ = constant_transport_coefficient::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<TransportCoefficientCommandBase> create_default_command(const Params& params) const override {
        if (params.physics.num_order_parameters < 0) {
            throw std::runtime_error("model transport constant requires nonnegative number of order parameters.");
        }

        auto cmd = std::make_shared<ConstantTransportCoefficientCommand>();
        cmd->type = constant_transport_coefficient::TYPE_NAME;
        cmd->num_order_parameters = params.physics.num_order_parameters;
        cmd->order_parameter_mobility.assign(static_cast<std::size_t>(cmd->num_order_parameters), 0.0);

        return cmd;
    }

    void update_command(TransportCoefficientCommandBase& command, const TransportCoefficientArgs& args, const Params& /*params*/) const override {
        auto* cmd = dynamic_cast<ConstantTransportCoefficientCommand*>(&command);
        if (!cmd) {
            throw std::runtime_error("ConstantTransportCoefficientStyle received incompatible command.");
        }

        for (const auto& kv : args.entries) {
            const std::string& key = kv.first;

            if (update_common_command(*cmd, key, kv.second)) {
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
            cmd->num_order_parameters,
            cmd->shear_viscosity,
            cmd->bulk_viscosity,
            cmd->order_parameter_mobility
        );
    }
};

#endif
