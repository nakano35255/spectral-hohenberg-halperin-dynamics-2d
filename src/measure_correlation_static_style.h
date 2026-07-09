#ifndef SHHD_MEASURE_CORRELATION_STATIC_STYLE_H
#define SHHD_MEASURE_CORRELATION_STATIC_STYLE_H

#include "measure_correlation_static.h"
#include "measure_registry.h"

#include <cctype>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace correlation_static_measure {
     inline const std::string TYPE_NAME = "correlation/static";
}

struct CorrelationStaticMeasureCommand : public MeasureCommandBase {
     int nevery = 10;
     int nblock = 10000;
     std::string file;
     CorrelationStaticAverageMode average_mode = CorrelationStaticAverageMode::Running;
     CorrelationStaticOutputMode output_mode = CorrelationStaticOutputMode::Shell;
     bool cross = false;
     std::vector<CorrelationStaticTarget> targets;

     void print(std::ostream& os) const override {
          const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
          os << "  " << std::left << std::setw(25)
               << label << ": id=" << id << " type=" << type;

          if (enabled) {
               os << " Nevery=" << nevery
                    << " Nblock=" << nblock
                    << " file=" << file
                    << " mode=" << (output_mode == CorrelationStaticOutputMode::TwoD ? "2d" : "shell")
                    << " average=" << (average_mode == CorrelationStaticAverageMode::Block ? "block" : "running")
                    << " cross=" << (cross ? "on" : "off")
                    << " target";

               for (const auto& target : targets) {
                    os << ' ' << target.name;
               }
          }

          os << "\n";
     }
};

class CorrelationStaticMeasureStyle : public MeasureStyle {
private:
     const std::string name_ = correlation_static_measure::TYPE_NAME;

     static bool parse_on_off(const std::string& value, const std::string& context) {
          if (value == "on") return true;
          if (value == "off") return false;
          throw std::runtime_error(context + " must be on|off.");
     }

     static CorrelationStaticAverageMode parse_average(const std::string& value) {
          if (value == "block") return CorrelationStaticAverageMode::Block;
          if (value == "running") return CorrelationStaticAverageMode::Running;
          throw std::runtime_error("correlation/static measure: average must be block|running.");
     }

     static CorrelationStaticOutputMode parse_mode(const std::string& value) {
          if (value == "2d") return CorrelationStaticOutputMode::TwoD;
          if (value == "shell") return CorrelationStaticOutputMode::Shell;
          throw std::runtime_error("correlation/static measure: mode must be 2d|shell.");
     }

     static bool parse_indexed_target(
          const std::string& value,
          const std::string& prefix,
          const std::string& suffix,
          int& component
     ) {
          if (value.rfind(prefix, 0) != 0) return false;
          if (value.size() <= prefix.size() + suffix.size()) return false;
          if (value.compare(value.size() - suffix.size(), suffix.size(), suffix) != 0) return false;

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

     static CorrelationStaticTarget parse_target(const std::string& target, const Params& params) {
          if (target == "rho") return {CorrelationStaticTargetKind::Rho, -1, target};
          if (target == "jx") return {CorrelationStaticTargetKind::Jx, -1, target};
          if (target == "jy") return {CorrelationStaticTargetKind::Jy, -1, target};

          int component = -1;
          if (parse_indexed_target(target, "psi[", "]", component)) {
               if (component < 0 || component >= params.physics.num_order_parameters) {
                    throw std::runtime_error("correlation/static target order-parameter index out of range: " + target);
               }
               return {CorrelationStaticTargetKind::Psi, component, target};
          }

          throw std::runtime_error("unknown correlation/static target: " + target);
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
          auto cmd = std::make_shared<CorrelationStaticMeasureCommand>();
          cmd->id = id;
          cmd->type = name_;
          cmd->enabled = enabled;

          if (!enabled) {
               return cmd;
          }

          cmd->nevery = std::stoi(args.get_required("nevery"));
          cmd->nblock = std::stoi(args.get_required("nblock"));
          cmd->file = args.get_required("file");
          cmd->output_mode = parse_mode(args.get_required("mode"));
          cmd->average_mode = parse_average(args.get_or("average", "running"));
          cmd->cross = parse_on_off(args.get_or("cross", "off"), "correlation/static cross");

          for (const std::string& target : args.targets) {
               cmd->targets.push_back(parse_target(target, params));
          }

          if (cmd->nevery <= 0) {
               throw std::runtime_error("correlation/static measure: nevery must be positive.");
          }
          if (cmd->nblock <= 0) {
               throw std::runtime_error("correlation/static measure: nblock must be positive.");
          }
          if (cmd->nblock % cmd->nevery != 0) {
               throw std::runtime_error("correlation/static measure: nblock must be divisible by nevery.");
          }
          if (cmd->file.empty()) {
               throw std::runtime_error("correlation/static measure: file is required.");
          }
          if (cmd->targets.empty()) {
               throw std::runtime_error("correlation/static measure requires target values.");
          }

          return cmd;
     }

     std::unique_ptr<Measure> create_measure(
          const Params& params,
          const Domain2D& domain,
          const Thermodynamics&,
          const FreeEnergy&,
          const TransportCoefficient&,
          std::shared_ptr<const MeasureCommandBase> command
     ) const override {
          auto static_cmd = std::dynamic_pointer_cast<const CorrelationStaticMeasureCommand>(command);
          if (!static_cmd) {
               throw std::runtime_error("CorrelationStaticMeasureStyle: invalid command type.");
          }

          return std::make_unique<CorrelationStaticMeasure>(params, domain, static_cmd);
     }
};

#endif