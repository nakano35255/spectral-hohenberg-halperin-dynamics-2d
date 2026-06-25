#ifndef SHHD_MEASURE_LONG_TIME_TAIL_STYLE_H
#define SHHD_MEASURE_LONG_TIME_TAIL_STYLE_H

#include "measure_long_time_tail.h"
#include "measure_registry.h"

#include <cctype>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace long_time_tail_measure {
     inline const std::string TYPE_NAME = "long_time_tail";
}

struct LongTimeTailMeasureCommand : public MeasureCommandBase {
     int nevery = 10;
     int nblock = 10000;
     std::string file;
     LongTimeTailAverageMode average_mode = LongTimeTailAverageMode::Running;
     bool cross = false;
     std::vector<LongTimeTailTarget> targets;

     void print(std::ostream& os) const override {
          const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
          os << "  " << std::left << std::setw(25)
               << label << ": id=" << id << " type=" << type;

          if (enabled) {
               os << " Nevery=" << nevery
                    << " Nblock=" << nblock
                    << " file=" << file
                    << " average=" << (average_mode == LongTimeTailAverageMode::Block ? "block" : "running")
                    << " cross=" << (cross ? "on" : "off")
                    << " target";

               for (const auto& target : targets) {
                    os << ' ' << target.name;
               }
          }

          os << "\n";
     }
};

class LongTimeTailMeasureStyle : public MeasureStyle {
private:
     const std::string name_ = long_time_tail_measure::TYPE_NAME;

     static bool parse_on_off(const std::string& value, const std::string& context) {
          if (value == "on") return true;
          if (value == "off") return false;
          throw std::runtime_error(context + " must be on|off.");
     }

     static LongTimeTailAverageMode parse_average(const std::string& value) {
          if (value == "block") return LongTimeTailAverageMode::Block;
          if (value == "running") return LongTimeTailAverageMode::Running;
          throw std::runtime_error("long_time_tail measure: average must be block|running.");
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

     static LongTimeTailTarget parse_target(const std::string& target, const Params& params) {
          if (target == "rho") return {LongTimeTailTargetKind::Rho, -1, target};
          if (target == "jx") return {LongTimeTailTargetKind::Jx, -1, target};
          if (target == "jy") return {LongTimeTailTargetKind::Jy, -1, target};

          int component = -1;
          if (parse_indexed_target(target, "psi[", "]", component)) {
               if (component < 0 || component >= params.physics.num_order_parameters) {
                    throw std::runtime_error("long_time_tail target order-parameter index out of range: " + target);
               }
               return {LongTimeTailTargetKind::Psi, component, target};
          }

          throw std::runtime_error("unknown long_time_tail target: " + target);
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
          auto cmd = std::make_shared<LongTimeTailMeasureCommand>();
          cmd->id = id;
          cmd->type = name_;
          cmd->enabled = enabled;

          if (!enabled) {
               return cmd;
          }

          cmd->nevery = std::stoi(args.get_required("nevery"));
          cmd->nblock = std::stoi(args.get_required("nblock"));
          cmd->file = args.get_required("file");
          cmd->average_mode = parse_average(args.get_or("average", "running"));
          cmd->cross = parse_on_off(args.get_or("cross", "off"), "long_time_tail cross");

          for (const std::string& target : args.targets) {
               cmd->targets.push_back(parse_target(target, params));
          }

          if (cmd->nevery <= 0) {
               throw std::runtime_error("long_time_tail measure: nevery must be positive.");
          }
          if (cmd->nblock <= 0) {
               throw std::runtime_error("long_time_tail measure: nblock must be positive.");
          }
          if (cmd->nblock % cmd->nevery != 0) {
               throw std::runtime_error("long_time_tail measure: nblock must be divisible by nevery.");
          }
          if (cmd->file.empty()) {
               throw std::runtime_error("long_time_tail measure: file is required.");
          }
          if (cmd->targets.empty()) {
               throw std::runtime_error("long_time_tail measure requires target values.");
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
          auto ltt_cmd = std::dynamic_pointer_cast<const LongTimeTailMeasureCommand>(command);
          if (!ltt_cmd) {
               throw std::runtime_error("LongTimeTailMeasureStyle: invalid command type.");
          }

          return std::make_unique<LongTimeTailMeasure>(params, domain, ltt_cmd);
     }
};

#endif