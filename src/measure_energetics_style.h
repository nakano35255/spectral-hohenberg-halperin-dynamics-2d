#ifndef SFI_MEASURE_ENERGETICS_STYLE_H
#define SFI_MEASURE_ENERGETICS_STYLE_H

#include "measure_energetics.h"
#include "measure_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace energetics_measure {
    inline const std::string TYPE_NAME = "energetics";
}

// ---------------------------------------------------------------------- //
struct EnergeticsMeasureCommand : public MeasureCommandBase {
    int nevery = 1000;
    std::string file;

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " Nevery=" << nevery
               << " file=" << file;
        }

        os << "\n";
    }
};
// ---------------------------------------------------------------------- //
class EnergeticsMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = energetics_measure::TYPE_NAME;

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MeasureCommandBase> parse_command(const std::string& id, bool enabled, const MeasureArgs& args, const Params& /*params*/) const override {
        auto cmd = std::make_shared<EnergeticsMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;

        if (!enabled) {
            return cmd;
        }

        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->file = args.get_required("file");

        if (cmd->nevery <= 0) {
            throw std::runtime_error("energetics measure: nevery must be positive.");
        }
        if (cmd->file.empty()) {
            throw std::runtime_error("energetics measure: file is required.");
        }

        return cmd;
    }

    std::unique_ptr<Measure> create_measure(
        const Params& params,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport,
        std::shared_ptr<const MeasureCommandBase> command
    ) const override {
        auto energetics_cmd = std::dynamic_pointer_cast<const EnergeticsMeasureCommand>(command);
        if (!energetics_cmd) {
            throw std::runtime_error("EnergeticsMeasureStyle: invalid command type.");
        }
        return std::make_unique<EnergeticsMeasure>(params, thermodynamics, free_energy, transport, energetics_cmd);
    }
};
// ---------------------------------------------------------------------- //

#endif
