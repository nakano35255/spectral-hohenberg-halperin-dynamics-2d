#ifndef SHHD_PACKAGE_PASSIVE_SCALAR_MEASURE_BUDGET_SPECTRUM_STYLE_H
#define SHHD_PACKAGE_PASSIVE_SCALAR_MEASURE_BUDGET_SPECTRUM_STYLE_H

#include "PASSIVE_SCALAR/measure_budget_spectrum.h"
#include "measure_registry.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace budget_spectrum_measure {
    inline const std::string TYPE_NAME = "budget/spectrum";
}

struct BudgetSpectrumMeasureCommand : public MeasureCommandBase {
    int component = 0;
    int nevery = 1000;
    int nblock = 1000;
    std::string file;
    BudgetSpectrumOutputMode output_mode = BudgetSpectrumOutputMode::TwoD;
    BudgetSpectrumAverageMode average_mode = BudgetSpectrumAverageMode::Block;

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " component=" << component
               << " Nevery=" << nevery
               << " Nblock=" << nblock
               << " file=" << file
               << " mode=" << (output_mode == BudgetSpectrumOutputMode::TwoD ? "2d" : "shell")
               << " average=" << (average_mode == BudgetSpectrumAverageMode::Block ? "block" : "running");
        }
        os << "\n";
    }
};

class BudgetSpectrumMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = budget_spectrum_measure::TYPE_NAME;

    static BudgetSpectrumOutputMode parse_mode(const std::string& value) {
        if (value == "2d") return BudgetSpectrumOutputMode::TwoD;
        if (value == "shell") return BudgetSpectrumOutputMode::Shell;
        throw std::runtime_error("budget/spectrum measure: mode must be 2d|shell.");
    }

    static BudgetSpectrumAverageMode parse_average(const std::string& value) {
        if (value == "block") return BudgetSpectrumAverageMode::Block;
        if (value == "running") return BudgetSpectrumAverageMode::Running;
        throw std::runtime_error("budget/spectrum measure: average must be block|running.");
    }

public:
    const std::string& type_name() const override { return name_; }

    std::shared_ptr<MeasureCommandBase> parse_command(
        const std::string& id,
        bool enabled,
        const MeasureArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<BudgetSpectrumMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;
        if (!enabled) return cmd;

        cmd->component = std::stoi(args.get_required("component"));
        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->nblock = std::stoi(args.get_required("nblock"));
        cmd->file = args.get_required("file");
        cmd->output_mode = parse_mode(args.get_required("mode"));
        cmd->average_mode = parse_average(args.get_required("average"));

        if (cmd->component < 0 || cmd->component >= params.physics.num_order_parameters) {
            throw std::runtime_error("budget/spectrum measure: component index out of range.");
        }
        if (cmd->nevery <= 0) throw std::runtime_error("budget/spectrum measure: nevery must be positive.");
        if (cmd->nblock <= 0) throw std::runtime_error("budget/spectrum measure: nblock must be positive.");
        if (cmd->nblock % cmd->nevery != 0) {
            throw std::runtime_error("budget/spectrum measure: nblock must be divisible by nevery.");
        }
        if (cmd->file.empty()) throw std::runtime_error("budget/spectrum measure: file is required.");
        return cmd;
    }

    std::unique_ptr<Measure> create_measure(
        const Params& params,
        const Domain2D& domain,
        const Thermodynamics&,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient,
        std::shared_ptr<const MeasureCommandBase> command
    ) const override {
        auto budget_cmd = std::dynamic_pointer_cast<const BudgetSpectrumMeasureCommand>(command);
        if (!budget_cmd) throw std::runtime_error("BudgetSpectrumMeasureStyle: invalid command type.");
        return std::make_unique<BudgetSpectrumMeasure>(
            params, domain, free_energy, transport_coefficient, budget_cmd
        );
    }
};

#endif