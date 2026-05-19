#ifndef MEASURE_SNAPSHOT_STYLE_H
#define MEASURE_SNAPSHOT_STYLE_H

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <iomanip>

#include "measure_registry.h"
#include "measure_snapshot.h"

namespace snapshot_measure {
    inline const std::string TYPE_NAME = "snapshot";
}

// ---------------------------------------------------------------------- //
struct SnapshotMeasureCommand : public MeasureCommandBase {
    int nevery = 100;
    std::string file;
    std::string space = "physical";

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " Nevery=" << nevery
               << " file=" << file
               << " space=" << space;
        }

        os << "\n";
    }
};
// ---------------------------------------------------------------------- //
class SnapshotMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = snapshot_measure::TYPE_NAME;

    static std::string canonical_space(const std::string& value) {
        if (value == "physical" || value == "real") {
            return "physical";
        }
        if (value == "spectral" || value == "fourier" || value == "spec") {
            return "spectral";
        }
        throw std::runtime_error("snapshot measure: space must be physical|spectral.");
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MeasureCommandBase> parse_command(const std::string& id, bool enabled, const MeasureArgs& args, const Params& /*params*/) const override {
        auto cmd = std::make_shared<SnapshotMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;

        if (!enabled) {
            return cmd;
        }

        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->file = args.get_required("file");
        cmd->space = canonical_space(args.get_or("space", "physical"));

        if (cmd->nevery <= 0) {
            throw std::runtime_error("snapshot measure: nevery must be positive.");
        }
        if (cmd->file.empty()) {
            throw std::runtime_error("snapshot measure: file is required.");
        }

        return cmd;
    }

    std::unique_ptr<Measure> create_measure(const Params& params, std::shared_ptr<const MeasureCommandBase> command) const override {
        auto snapshot_cmd = std::dynamic_pointer_cast<const SnapshotMeasureCommand>(command);
        if (!snapshot_cmd) {
            throw std::runtime_error("SnapshotMeasureStyle: invalid command type.");
        }
        return std::make_unique<SnapshotMeasure>(params, snapshot_cmd);
    }
};
// ---------------------------------------------------------------------- //
#endif
