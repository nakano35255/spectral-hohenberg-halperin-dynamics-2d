#ifndef SHHD_MEASURE_YOKOTA_GREEN_KUBO_STYLE_H
#define SHHD_MEASURE_YOKOTA_GREEN_KUBO_STYLE_H

#include "measure_registry.h"
#include "measure_yokota_green_kubo.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace yokota_green_kubo_measure {
    inline const std::string TYPE_NAME = "yokota_green_kubo";
}

// ---------------------------------------------------------------------- //
struct YokotaGreenKuboMeasureCommand : public MeasureCommandBase {
    int nevery = 1000;
    int nblock = 1000;
    std::string file;
    std::string mode = "diagonal";

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " Nevery=" << nevery
               << " Nblock=" << nblock
               << " file=" << file
               << " mode=" << mode;
        }

        os << "\n";
    }
};
// ---------------------------------------------------------------------- //
class YokotaGreenKuboMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = yokota_green_kubo_measure::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MeasureCommandBase> parse_command(
        const std::string& id,
        bool enabled,
        const MeasureArgs& args,
        const Params& /*params*/
    ) const override {
        auto cmd = std::make_shared<YokotaGreenKuboMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;

        if (!enabled) {
            return cmd;
        }

        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->nblock = std::stoi(args.get_required("nblock"));
        cmd->file = args.get_required("file");
        cmd->mode = args.get_or("mode", "diagonal");

        if (cmd->nevery <= 0) {
            throw std::runtime_error("yokota_green_kubo measure: nevery must be positive.");
        }
        if (cmd->nblock <= 0) {
            throw std::runtime_error("yokota_green_kubo measure: nblock must be positive.");
        }
        if (cmd->nblock % cmd->nevery != 0) {
            throw std::runtime_error("yokota_green_kubo measure: nblock must be divisible by nevery.");
        }
        if (cmd->file.empty()) {
            throw std::runtime_error("yokota_green_kubo measure: file is required.");
        }
        if (cmd->mode != "diagonal" && cmd->mode != "all") {
            throw std::runtime_error("yokota_green_kubo measure: mode must be diagonal|all.");
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
        auto ykgk_cmd = std::dynamic_pointer_cast<const YokotaGreenKuboMeasureCommand>(command);
        if (!ykgk_cmd) {
            throw std::runtime_error("YokotaGreenKuboMeasureStyle: invalid command type.");
        }
        return std::make_unique<YokotaGreenKuboMeasure>(params, domain, ykgk_cmd);
    }
};
// ---------------------------------------------------------------------- //

#endif
