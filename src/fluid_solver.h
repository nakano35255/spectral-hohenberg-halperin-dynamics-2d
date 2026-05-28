#ifndef SFI_FLUID_SOLVER_H
#define SFI_FLUID_SOLVER_H

#include "domain.h"
#include "simulationinfo.h"
#include "state.h"
#include "spectral_mask.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "fcalculator.h"
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

     const Thermodynamics& thermodynamics_;
     const FreeEnergy& free_energy_;
     const TransportCoefficient& transport_coefficient_;

     SpectralMask2D spectral_mask_;
     FCalculator fcalculator_;
     std::unique_ptr<TimeIntegrator> time_integrator_;
     RHSOperators rhs_;

     static std::unique_ptr<TimeIntegrator> create_time_integrator(const Params& params, const Domain2D& domain, const SpectralMask2D& spectral_mask) {
          const std::string& type = params.runtime.time_evolution_type;

          if (type == "euler/compressible") return std::make_unique<EulerCompressible>(domain, params, spectral_mask);
          if (type == "euler/quiescent") return std::make_unique<EulerQuiescent>(domain, params, spectral_mask);
          if (type == "euler/incompressible") return std::make_unique<EulerIncompressible>(domain, params, spectral_mask);
          if (type == "srk3/compressible") return std::make_unique<SRK3Compressible>(domain, params, spectral_mask);
          if (type == "srk3/quiescent") return std::make_unique<SRK3Quiescent>(domain, params, spectral_mask);
          if (type == "srk3/incompressible") return std::make_unique<SRK3Incompressible>(domain, params, spectral_mask);

          throw std::runtime_error("unknown time_evolution type: " + type);
     }


public:
     FluidSolver(
          const Params& params,
          const Domain2D& domain,
          const Thermodynamics& thermodynamics,
          const FreeEnergy& free_energy,
          const TransportCoefficient& transport_coefficient
     )
          : params_(params),
          domain_(domain),
          thermodynamics_(thermodynamics),
          free_energy_(free_energy),
          transport_coefficient_(transport_coefficient),
          spectral_mask_(params_, domain_),
          fcalculator_(params_, domain_, spectral_mask_, thermodynamics_, free_energy_, transport_coefficient_),
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
          if (params_.fix.noise.enabled) {
               if (params_.physics.num_order_parameters > 0) {
                    rhs_.psi_sto = [this](int order_parameter, const State& state, Complex* out) {
                         fcalculator_.psi_sto(order_parameter, state, out);
                    };
               }

               const bool quiescent = params_.runtime.time_evolution_type == "euler/quiescent" || params_.runtime.time_evolution_type == "srk3/quiescent";
               if (!quiescent) {
                    rhs_.j_sto = [this](const State& state, Complex* out_jx, Complex* out_jy) {
                         fcalculator_.j_sto(state, out_jx, out_jy);
                    };
               }
          }
     }

     void step(State& state, double time) {
          time_integrator_->step(state, time, rhs_);
     }
};

#endif
