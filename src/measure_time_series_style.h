#ifndef SHHD_MEASURE_TIME_SERIES_STYLE_H
#define SHHD_MEASURE_TIME_SERIES_STYLE_H

#include "measure_registry.h"
#include "measure_time_series.h"

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cctype>

namespace time_series_measure {
    inline const std::string TYPE_NAME = "time_series";
}

struct TimeSeriesMeasureCommand : public MeasureCommandBase {
    int nevery = 1000;
    std::string file;
    std::vector<std::string> targets;
    std::vector<TimeSeriesColumn> columns;

    void print(std::ostream& os) const override {
        const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
        os << "  " << std::left << std::setw(25)
           << label << ": id=" << id << " type=" << type;
        if (enabled) {
            os << " Nevery=" << nevery
               << " file=" << file;
            os << " target";
            for (const TimeSeriesColumn& column : columns) {
                os << ' ' << column.name;
            }
        }
        os << "\n";
    }
};

class TimeSeriesMeasureStyle : public MeasureStyle {
private:
    const std::string name_ = time_series_measure::TYPE_NAME;

    static bool parse_indexed_suffix(const std::string& value, const std::string& prefix, const std::string& suffix, int& component) {
        if (value.rfind(prefix, 0) != 0) {
            return false;
        }
        if (value.size() <= prefix.size() + suffix.size()) {
            return false;
        }
        if (value.compare(value.size() - suffix.size(), suffix.size(), suffix) != 0) {
            return false;
        }

        const std::size_t begin = prefix.size();
        const std::size_t end = value.size() - suffix.size();
        for (std::size_t i = begin; i < end; ++i) {
            if (!std::isdigit(static_cast<unsigned char>(value[i]))) {
                return false;
            }
        }

        component = std::stoi(value.substr(begin, end - begin));
        return true;
    }

    static TimeSeriesColumn parse_column(const std::string& target, const Params& params) {
        if (target == "E_T") return {TimeSeriesColumnKind::EnergyTotal, -1, target};
        if (target == "E_K") return {TimeSeriesColumnKind::EnergyKinetic, -1, target};
        if (target == "E_psi") return {TimeSeriesColumnKind::EnergyFree, -1, target};
        if (target == "E_C") return {TimeSeriesColumnKind::EnergyCompressibility, -1, target};
        if (target == "pi_xx") return {TimeSeriesColumnKind::MomentumFluxXX, -1, target};
        if (target == "pi_xy") return {TimeSeriesColumnKind::MomentumFluxXY, -1, target};
        if (target == "pi_yy") return {TimeSeriesColumnKind::MomentumFluxYY, -1, target};
        if (target == "|rho|^2") return {TimeSeriesColumnKind::DensityL2, -1, target};
        if (target == "|d_rho|^2") return {TimeSeriesColumnKind::DensityGradientL2, -1, target};
        if (target == "|j|^2") return {TimeSeriesColumnKind::MomentumL2, -1, target};
        if (target == "|d_j|^2") return {TimeSeriesColumnKind::MomentumGradientL2, -1, target};

        int component = -1;
        if (parse_indexed_suffix(target, "Jpsi[", "]_x", component)) {
            if (component < 0 || component >= params.physics.num_order_parameters) {
                throw std::runtime_error("time_series target order-parameter index out of range: " + target);
            }
            return {TimeSeriesColumnKind::OrderParameterFluxX, component, target};
        }

        if (parse_indexed_suffix(target, "Jpsi[", "]_y", component)) {
            if (component < 0 || component >= params.physics.num_order_parameters) {
                throw std::runtime_error("time_series target order-parameter index out of range: " + target);
            }
            return {TimeSeriesColumnKind::OrderParameterFluxY, component, target};
        }

        if (parse_indexed_suffix(target, "|psi[", "]|^2", component)) {
            if (component < 0 || component >= params.physics.num_order_parameters) {
                throw std::runtime_error("time_series target order-parameter index out of range: " + target);
            }
            return {TimeSeriesColumnKind::OrderParameterL2, component, target};
        }

        if (parse_indexed_suffix(target, "|d_psi[", "]|^2", component)) {
            if (component < 0 || component >= params.physics.num_order_parameters) {
                throw std::runtime_error("time_series target order-parameter index out of range: " + target);
            }
            return {TimeSeriesColumnKind::OrderParameterGradientL2, component, target};
        }

        throw std::runtime_error("unknown time_series target: " + target);
    }

public:
    const std::string& type_name() const override {
        return name_;
    }

    std::shared_ptr<MeasureCommandBase> parse_command(
        const std::string& id,
        bool enabled,
        const MeasureArgs& args,
        const Params& params
    ) const override {
        auto cmd = std::make_shared<TimeSeriesMeasureCommand>();
        cmd->id = id;
        cmd->type = name_;
        cmd->enabled = enabled;

        if (!enabled) {
            return cmd;
        }

        cmd->nevery = std::stoi(args.get_required("nevery"));
        cmd->file = args.get_required("file");
        cmd->targets = args.targets;

        for (const std::string& target : args.targets) {
            cmd->columns.push_back(parse_column(target, params));
        }

        if (cmd->targets.empty()) {
            throw std::runtime_error("time_series measure requires target values.");
        }
        if (cmd->nevery <= 0) {
            throw std::runtime_error("time_series measure: nevery must be positive.");
        }
        if (cmd->file.empty()) {
            throw std::runtime_error("time_series measure: file is required.");
        }
        if (cmd->columns.empty()) {
            throw std::runtime_error("time_series measure requires target values.");
        }
        return cmd;
    }

    std::unique_ptr<Measure> create_measure(
        const Params& params,
        const Domain2D& domain,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient,
        std::shared_ptr<const MeasureCommandBase> command
    ) const override {
        auto time_series_cmd = std::dynamic_pointer_cast<const TimeSeriesMeasureCommand>(command);
        if (!time_series_cmd) {
            throw std::runtime_error("TimeSeriesMeasureStyle: invalid command type.");
        }

        return std::make_unique<TimeSeriesMeasure>(
            params, domain, thermodynamics, free_energy, transport_coefficient, time_series_cmd
        );
    }
};

#endif