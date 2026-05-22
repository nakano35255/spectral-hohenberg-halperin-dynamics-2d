#ifndef SFI_FLUID_SOLVER_H
#define SFI_FLUID_SOLVER_H

#include "domain.h"
#include "fcalculator.h"
#include "model_thermodynamics.h"
#include "model_thermodynamics_registry_builtin.h"
#include "model_transport_coefficient.h"
#include "model_transport_coefficient_registry_builtin.h"
#include "simulationinfo.h"
#include "state.h"
#include "time_integrator.h"
#include "time_integrator_euler_compressible.h"
#include "time_integrator_euler_incompressible.h"
#include "time_integrator_euler_quiescent.h"
#include "time_integrator_srk3_compressible.h"
#include "time_integrator_srk3_incompressible.h"
#include "time_integrator_srk3_quiescent.h"

#include <memory>
#include <string>

class FluidSolver {
private:
     const Params& params_;
     const Domain2D& domain_;

     ThermodynamicsRegistry thermodynamics_registry_;
     TransportCoefficientRegistry transport_coefficient_registry_;

     std::unique_ptr<Thermodynamics> thermodynamics_;
     std::unique_ptr<TransportCoefficient> transport_coefficient_;
     FCalculator fcalculator_;
     std::unique_ptr<TimeIntegrator> time_integrator_;
     RHSOperators rhs_;

     static bool is_incompressible_time_evolution(const std::string& type) {
          return type == "euler/incompressible" || type == "srk3/incompressible";
     }

     static std::unique_ptr<Thermodynamics> create_thermodynamics(
          const Params& params,
          const ThermodynamicsRegistry& registry
     ) {
          const std::string& thermo_type = params.physics.thermo->type;
          const ThermodynamicsStyle& style = registry.get_thermo(thermo_type);
          return style.create(params, params.physics.thermo);
     }

     static std::unique_ptr<TransportCoefficient> create_transport_coefficient(
          const Params& params,
          const TransportCoefficientRegistry& registry
     ) {
          const std::string& transport_type = params.physics.transport->type;
          const TransportCoefficientStyle& style = registry.get_transport(transport_type);
          return style.create(params, params.physics.transport);
     }

     static std::unique_ptr<TimeIntegrator> create_time_integrator(
          const Params& params,
          const Domain2D& domain
     ) {
          const std::string& type = params.runtime.time_evolution_type;

          if (type == "euler/compressible") {
               return std::make_unique<EulerCompressible>(domain, params);
          }

          if (type == "euler/quiescent") {
               return std::make_unique<EulerQuiescent>(domain, params);
          }

          if (type == "euler/incompressible") {
               return std::make_unique<EulerIncompressible>(domain, params);
          }

          if (type == "srk3/compressible") {
               return std::make_unique<SRK3Compressible>(domain, params);
          }

          if (type == "srk3/quiescent") {
               return std::make_unique<SRK3Quiescent>(domain, params);
          }

          if (type == "srk3/incompressible") {
               return std::make_unique<SRK3Incompressible>(domain, params);
          }

          throw std::runtime_error("unknown time_evolution type: " + type);
     }


public:
     FluidSolver(const Params& params, const Domain2D& domain)
          : params_(params),
            domain_(domain),
            thermodynamics_registry_(build_thermodynamics_registry()),
            transport_coefficient_registry_(build_transport_coefficient_registry()),
            thermodynamics_(create_thermodynamics(params_, thermodynamics_registry_)),
            transport_coefficient_(create_transport_coefficient(params_, transport_coefficient_registry_)),
            fcalculator_(params_, domain_, *thermodynamics_, *transport_coefficient_),
            time_integrator_(create_time_integrator(params_, domain_)) {
          rhs_.density_det = [this](int component, const State& state, Complex* out, double time) {
               fcalculator_.density_det(component, state, out, time);
          };
          if (is_incompressible_time_evolution(params_.runtime.time_evolution_type)) {
               rhs_.momentum_det = [this](const State& state, Complex* out_jx, Complex* out_jy, double time) {
                    fcalculator_.momentum_det_incompressible(state, out_jx, out_jy, time);
               };
          } else {
               rhs_.momentum_det = [this](const State& state, Complex* out_jx, Complex* out_jy, double time) {
                    fcalculator_.momentum_det(state, out_jx, out_jy, time);
               };
          }
     }

     void step(State& state, double time) {
          time_integrator_->step(state, time, rhs_);
     }
};

#endif
