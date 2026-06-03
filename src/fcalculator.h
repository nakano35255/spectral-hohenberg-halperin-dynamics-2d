#ifndef SHHD_FCALCULATOR_H
#define SHHD_FCALCULATOR_H

#include "domain.h"
#include "fcalculator_dynamics_mode.h"
#include "fcalculator_field_requests_and_layout.h"
#include "fcalculator_workspace.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "simulationinfo.h"
#include "spectral_mask.h"
#include "state.h"
#include "buffer_flux.h"

#include <cstddef>
#include <random>
#include <vector>

class FCalculator {
private:
    const Params& params_;
    const Domain2D& domain_;
    const SpectralMask2D& spectral_mask_;
    const Thermodynamics& thermodynamics_;
    const FreeEnergy& free_energy_;
    const TransportCoefficient& transport_coefficient_;

    DynamicsMode dynamics_mode_;
    PhysicalRHSFieldRequests physical_field_requests_;

    int num_order_parameters_ = 0;
    std::size_t local_spectral_size_ = 0;

    mutable FCalculatorWorkspace workspace_;

    mutable std::mt19937_64 noise_rng_;
    mutable std::normal_distribution<double> normal_noise_;

    static int temporary_field_capacity(const Params& params);
    void clear(Complex* out) const;

    void add_order_parameter_linear_term(int order_parameter, double mobility, const State& current, Complex* out, FluxBuffer* flux) const;
    void add_order_parameter_physical_terms(int order_parameter, double mobility, const State& current, Complex* out, double time, FluxBuffer* flux) const;
    void add_linear_pressure_term(double pressure_coefficient, const State& current, Complex* out_jx, Complex* out_jy, FluxBuffer* flux) const;
    void add_incompressible_viscous_term(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy, FluxBuffer* flux) const;
    void add_compressible_viscous_term(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy, double time, FluxBuffer* flux) const;
    void add_momentum_physical_terms(double eta, double zeta, const State& current, Complex* out_jx, Complex* out_jy, double time, FluxBuffer* flux) const;

    void add_momentum_sine_force_terms(Complex* out_jx, Complex* out_jy) const;
    void add_order_parameter_sine_force_terms(int order_parameter, Complex* out) const;
    void add_order_parameter_gradient_force_terms(int order_parameter, const State& current, Complex* out, double time) const;
    
public:
    FCalculator(
        const Params& params,
        const Domain2D& domain,
        const SpectralMask2D& spectral_mask,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient
    );

    void rho_det(const State& current, Complex* out, double t, FluxBuffer* flux) const;
    void psi_det(int order_parameter, const State& current, Complex* out, double t, FluxBuffer* flux) const;
    void j_det(const State& current, Complex* out_jx, Complex* out_jy, double t, FluxBuffer* flux) const;

    void psi_sto(int order_parameter, const State& current, Complex* out, FluxBuffer* flux) const;
    void j_sto(const State& current, Complex* out_jx, Complex* out_jy, FluxBuffer* flux) const;
};

#endif
