#ifndef SFI_FLUID_SOLVER_H
#define SFI_FLUID_SOLVER_H

#include "domain.h"
#include "fcalculator.h"
#include "model_free_energy.h"
#include "model_free_energy_null.h"
#include "model_free_energy_registry_builtin.h"
#include "model_thermodynamics.h"
#include "model_thermodynamics_null.h"
#include "model_thermodynamics_registry_builtin.h"
#include "model_transport_coefficient.h"
#include "model_transport_coefficient_registry_builtin.h"
#include "simulationinfo.h"
#include "state.h"
#include "spectral_mask.h"
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
     FreeEnergyRegistry free_energy_registry_;
     TransportCoefficientRegistry transport_coefficient_registry_;

     std::unique_ptr<Thermodynamics> thermodynamics_;
     std::unique_ptr<FreeEnergy> free_energy_;
     std::unique_ptr<TransportCoefficient> transport_coefficient_;

     SpectralMask2D spectral_mask_;
     FCalculator fcalculator_;
     std::unique_ptr<TimeIntegrator> time_integrator_;
     RHSOperators rhs_;

     static std::unique_ptr<Thermodynamics> create_thermodynamics(
          const Params& params,
          const ThermodynamicsRegistry& registry
     ) {
          if (!params.physics.thermo) {
               return std::make_unique<NullThermodynamics>();
          }

          const std::string& thermo_type = params.physics.thermo->type;
          const ThermodynamicsStyle& style = registry.get_thermo(thermo_type);
          return style.create(params, params.physics.thermo);
     }

     static std::unique_ptr<FreeEnergy> create_free_energy(
          const Params& params,
          const FreeEnergyRegistry& registry
     ) {
          if (!params.physics.free_energy) {
               return std::make_unique<NullFreeEnergy>(params.physics.num_order_parameters);
          }

          const std::string& free_energy_type = params.physics.free_energy->type;
          const FreeEnergyStyle& style = registry.get_free_energy(free_energy_type);
          return style.create(params, params.physics.free_energy);
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
          const Domain2D& domain,
          const SpectralMask2D& spectral_mask
     ) {
          const std::string& type = params.runtime.time_evolution_type;

          if (type == "euler/compressible") {
               return std::make_unique<EulerCompressible>(domain, params, spectral_mask);
          }

          if (type == "euler/quiescent") {
               return std::make_unique<EulerQuiescent>(domain, params, spectral_mask);
          }

          if (type == "euler/incompressible") {
               return std::make_unique<EulerIncompressible>(domain, params, spectral_mask);
          }

          if (type == "srk3/compressible") {
               return std::make_unique<SRK3Compressible>(domain, params, spectral_mask);
          }

          if (type == "srk3/quiescent") {
               return std::make_unique<SRK3Quiescent>(domain, params, spectral_mask);
          }

          if (type == "srk3/incompressible") {
               return std::make_unique<SRK3Incompressible>(domain, params, spectral_mask);
          }

          throw std::runtime_error("unknown time_evolution type: " + type);
     }


public:
     FluidSolver(const Params& params, const Domain2D& domain)
          : params_(params),
            domain_(domain),
            thermodynamics_registry_(build_thermodynamics_registry()),
            free_energy_registry_(build_free_energy_registry()),
            transport_coefficient_registry_(build_transport_coefficient_registry()),
            thermodynamics_(create_thermodynamics(params_, thermodynamics_registry_)),
            free_energy_(create_free_energy(params_, free_energy_registry_)),
            transport_coefficient_(create_transport_coefficient(params_, transport_coefficient_registry_)),
            spectral_mask_(domain_, params_),
            fcalculator_(params_, domain_, spectral_mask_, *thermodynamics_, *free_energy_, *transport_coefficient_),
            time_integrator_(create_time_integrator(params_, domain_, spectral_mask_))
     {
          rhs_.rho_det = [this](const State& state, Complex* out, double time) {
               fcalculator_.rho_det(state, out, time);
          };
          rhs_.psi_det = [this](int order_parameter, const State& state, Complex* out, double time) {
               fcalculator_.psi_det(order_parameter, state, out, time);
          };
          rhs_.j_det = [this](const State& state, Complex* out_jx, Complex* out_jy, double time) {
               fcalculator_.j_det(state, out_jx, out_jy, time);
          };

          if (fix_enabled(params_.fix.flags, FixFlag::Noise)) {
               rhs_.psi_sto = [this](int order_parameter, const State& state, Complex* out) {
                    fcalculator_.psi_sto(order_parameter, state, out);
               };
          }
     }

     void step(State& state, double time) {
          time_integrator_->step(state, time, rhs_);
     }
};

#endif
