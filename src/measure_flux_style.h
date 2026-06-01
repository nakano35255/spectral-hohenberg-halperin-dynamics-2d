#ifndef SHHD_MEASURE_FLUX_STYLE_H
#define SHHD_MEASURE_FLUX_STYLE_H

#include "measure_flux.h"
#include "measure_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace flux_measure {
    inline const std::string TYPE_NAME = "flux";
}

// ---------------------------------------------------------------------- //
struct FluxMeasureCommand : public MeasureCommandBase {
    int nevery = 1000;
    std::string file;
    FluxRequest request;
    std::string target = "all";

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " Nevery=" << nevery
               << " file=" << file
               << " target=" << target;
        }

        os << "\n";
    }
};
// ---------------------------------------------------------------------- //
class FluxMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = flux_measure::TYPE_NAME;

    static FluxRequest parse_target(const std::string& target, const Params& params) {
        FluxRequest request;

        if (target == "all") {
            request.density = true;
            request.order_parameter = params.physics.num_order_parameters > 0;
            request.momentum = true;
            return request;
        }
        if (target == "density" || target == "rho") {
            request.density = true;
            return request;
        }
        if (target == "order_parameter" || target == "psi") {
            if (params.physics.num_order_parameters <= 0) {
                throw std::runtime_error("flux measure: target order_parameter requires order_parameters > 0.");
            }
            request.order_parameter = true;
            return request;
        }
        if (target == "momentum" || target == "j") {
            request.momentum = true;
            return request;
        }

        throw std::runtime_error("flux measure: target must be all|density|order_parameter|momentum.");
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MeasureCommandBase> parse_command(const std::string& id, bool enabled, const MeasureArgs& args, const Params& params) const override {
        auto cmd = std::make_shared<FluxMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;

        if (!enabled) {
            return cmd;
        }

        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->file = args.get_required("file");
        cmd->target = args.get_or("target", "all");
        cmd->request = parse_target(cmd->target, params);

        if (cmd->nevery <= 0) {
            throw std::runtime_error("flux measure: nevery must be positive.");
        }
        if (cmd->file.empty()) {
            throw std::runtime_error("flux measure: file is required.");
        }

        return cmd;
    }

    std::unique_ptr<Measure> create_measure(
        const Params& params,
        const Domain2D& domain,
        const Thermodynamics& /*thermodynamics*/,
        const FreeEnergy& /*free_energy*/,
        const TransportCoefficient& /*transport_coefficient*/,
        std::shared_ptr<const MeasureCommandBase> command
    ) const override {
        auto flux_cmd = std::dynamic_pointer_cast<const FluxMeasureCommand>(command);
        if (!flux_cmd) {
            throw std::runtime_error("FluxMeasureStyle: invalid command type.");
        }
        return std::make_unique<FluxMeasure>(params, domain, flux_cmd);
    }
};
// ---------------------------------------------------------------------- //

#endif
