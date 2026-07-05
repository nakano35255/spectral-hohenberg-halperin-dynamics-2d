#ifndef SHHD_IMPLICIT_OPERATOR_H
#define SHHD_IMPLICIT_OPERATOR_H

#include "domain.h"
#include "fcalculator_dynamics_mode.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "simulationinfo.h"
#include "spectral_mask.h"
#include "state.h"

#include <algorithm>
#include <cstddef>

class ImplicitOperator {
private:
     const Params& params_;
     const Domain2D& domain_;
     const SpectralMask2D& spectral_mask_;
     const Thermodynamics& thermodynamics_;
     const FreeEnergy& free_energy_;
     const TransportCoefficient& transport_coefficient_;

     DynamicsMode dynamics_mode_;

     mutable bool reference_density_ready_ = false;
     mutable double reference_density_ = 0.0;
     double reference_density(const State& state) const;

     void solve_order_parameter_fields(State& dst, const State& src, double alpha) const;
     void solve_incompressible_hydrodynamic_fields(State& dst, const State& src, double alpha) const;
     void solve_compressible_hydrodynamic_fields(State& dst, const State& src, double alpha) const;

public:
     ImplicitOperator(
          const Params& params,
          const Domain2D& domain,
          const SpectralMask2D& spectral_mask,
          const Thermodynamics& thermodynamics,
          const FreeEnergy& free_energy,
          const TransportCoefficient& transport_coefficient
     )
          : params_(params),
               domain_(domain),
               spectral_mask_(spectral_mask),
               thermodynamics_(thermodynamics),
               free_energy_(free_energy),
               transport_coefficient_(transport_coefficient),
               dynamics_mode_(parse_dynamics_mode(params.runtime.time_evolution_type)) {}

     void apply_inverse(State& dst, const State& src, double alpha) const;

};

#endif