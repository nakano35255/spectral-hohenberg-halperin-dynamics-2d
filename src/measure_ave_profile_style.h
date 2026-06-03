#ifndef SHHD_MEASURE_AVE_PROFILE_STYLE_H
#define SHHD_MEASURE_AVE_PROFILE_STYLE_H

#include "measure_ave_profile.h"
#include "measure_registry.h"

#include <cctype>
#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace ave_profile_measure {
     inline const std::string TYPE_NAME = "ave/profile";
}

struct AveProfileMeasureCommand : public MeasureCommandBase {
     int axis = 1;
     int nevery = 1000;
     int nblock = 1000;
     std::string file;
     AveProfileAverageMode average_mode = AveProfileAverageMode::Block;
     std::vector<AveProfileTarget> targets;

     void print(std::ostream& os) const override {
          const std::string label = std::string("Measure ") + (enabled ? "on" : "off");
          os << "  " << std::left << std::setw(25)
               << label << ": id=" << id << " type=" << type;
          if (enabled) {
               os << " axis=" << (axis == 0 ? "x" : "y")
                    << " Nevery=" << nevery
                    << " Nblock=" << nblock
                    << " file=" << file
                    << " average=" << (average_mode == AveProfileAverageMode::Block ? "block" : "running")
                    << " target";
               for (const auto& target : targets) {
                    os << ' ' << target.name;
               }
          }
          os << "\n";
     }
};

class AveProfileMeasureStyle : public MeasureStyle {
private:
     const std::string name_ = ave_profile_measure::TYPE_NAME;

     static int parse_axis(const std::string& value) {
          if (value == "x" || value == "0") return 0;
          if (value == "y" || value == "1") return 1;
          throw std::runtime_error("ave/profile measure: axis must be x|y|0|1.");
     }

     static AveProfileAverageMode parse_average(const std::string& value) {
          if (value == "block") return AveProfileAverageMode::Block;
          if (value == "running") return AveProfileAverageMode::Running;
          throw std::runtime_error("ave/profile measure: average must be block|running.");
     }

     static bool parse_indexed_suffix(const std::string& value, const std::string& prefix, const std::string& suffix, int& component) {
          if (value.rfind(prefix, 0) != 0) return false;
          if (value.size() <= prefix.size() + suffix.size()) return false;
          if (value.compare(value.size() - suffix.size(), suffix.size(), suffix) != 0) return false;

          const std::size_t begin = prefix.size();
          const std::size_t end = value.size() - suffix.size();
          for (std::size_t i = begin; i < end; ++i) {
               if (!std::isdigit(static_cast<unsigned char>(value[i]))) return false;
          }

          component = std::stoi(value.substr(begin, end - begin));
          return true;
     }

     static AveProfileTarget parse_target(const std::string& target, const Params& params) {
          if (target == "rho") return {AveProfileTargetKind::Rho, -1, target};
          if (target == "jx") return {AveProfileTargetKind::Jx, -1, target};
          if (target == "jy") return {AveProfileTargetKind::Jy, -1, target};
          if (target == "vx") return {AveProfileTargetKind::Vx, -1, target};
          if (target == "vy") return {AveProfileTargetKind::Vy, -1, target};

          if (target == "pi_xx" || target == "pixx") return {AveProfileTargetKind::MomentumFluxXX, -1, "pi_xx"};
          if (target == "pi_xy" || target == "pixy") return {AveProfileTargetKind::MomentumFluxXY, -1, "pi_xy"};
          if (target == "pi_yy" || target == "piyy") return {AveProfileTargetKind::MomentumFluxYY, -1, "pi_yy"};

          int component = -1;
          if (parse_indexed_suffix(target, "psi[", "]", component)) {
               if (component < 0 || component >= params.physics.num_order_parameters) {
                    throw std::runtime_error("ave/profile target order-parameter index out of range: " + target);
               }
               return {AveProfileTargetKind::Psi, component, target};
          }
          if (parse_indexed_suffix(target, "Jpsi[", "]_x", component)) {
               if (component < 0 || component >= params.physics.num_order_parameters) {
                    throw std::runtime_error("ave/profile target order-parameter index out of range: " + target);
               }
               return {AveProfileTargetKind::OrderParameterFluxX, component, target};
          }
          if (parse_indexed_suffix(target, "Jpsi[", "]_y", component)) {
               if (component < 0 || component >= params.physics.num_order_parameters) {
                    throw std::runtime_error("ave/profile target order-parameter index out of range: " + target);
               }
               return {AveProfileTargetKind::OrderParameterFluxY, component, target};
          }

          throw std::runtime_error("unknown ave/profile target: " + target);
     }

public:
     const std::string& type_name() const override { return name_; }

     std::shared_ptr<MeasureCommandBase> parse_command(
          const std::string& id,
          bool enabled,
          const MeasureArgs& args,
          const Params& params
     ) const override {
          auto cmd = std::make_shared<AveProfileMeasureCommand>();
          cmd->id = id;
          cmd->type = name_;
          cmd->enabled = enabled;

          if (!enabled) return cmd;

          cmd->axis = parse_axis(args.get_required("axis"));
          cmd->nevery = std::stoi(args.get_required("nevery"));
          cmd->nblock = std::stoi(args.get_required("nblock"));
          cmd->file = args.get_required("file");
          cmd->average_mode = parse_average(args.get_required("average"));
          for (const std::string& target : args.targets) {
               cmd->targets.push_back(parse_target(target, params));
          }

          if (cmd->nevery <= 0) {
               throw std::runtime_error("ave/profile measure: nevery must be positive.");
          }
          if (cmd->nblock <= 0) {
               throw std::runtime_error("ave/profile measure: nblock must be positive.");
          }
          if (cmd->nblock % cmd->nevery != 0) {
               throw std::runtime_error("ave/profile measure: nblock must be divisible by nevery.");
          }
          if (cmd->file.empty()) {
               throw std::runtime_error("ave/profile measure: file is required.");
          }
          if (cmd->targets.empty()) {
               throw std::runtime_error("ave/profile measure requires target values.");
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
          auto profile_cmd = std::dynamic_pointer_cast<const AveProfileMeasureCommand>(command);
          if (!profile_cmd) throw std::runtime_error("AveProfileMeasureStyle: invalid command type.");
          return std::make_unique<AveProfileMeasure>(params, domain, profile_cmd);
     }
};

#endif