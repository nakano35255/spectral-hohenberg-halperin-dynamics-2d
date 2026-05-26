#ifndef SFI_FCALCULATOR_H
#define SFI_FCALCULATOR_H

#include "domain.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "simulationinfo.h"
#include "state.h"
#include "spectral_mask.h"
#include "fcalculator_workspace.h"

#include <cstddef>
#include <cstdint>
#include <vector>
#include <random>

class FCalculator {
private:
    const Params& params_;
    const Domain2D& domain_;
    const SpectralMask2D& spectral_mask_;
    const Thermodynamics& thermodynamics_;
    const FreeEnergy& free_energy_;
    const TransportCoefficient& transport_coefficient_;

    int num_order_parameters_ = 0;
    std::size_t local_spectral_size_ = 0;

    // for noise
    mutable std::mt19937_64 noise_rng_;
    mutable std::normal_distribution<double> normal_noise_;

    // for ifft plan
    bool is_quiescent_time_evolution() const;
    bool is_incompressible_time_evolution() const;
    bool is_compressible_time_evolution() const;
    PhysicalRHSPlan physical_rhs_plan_;
    PhysicalRHSPlan make_physical_rhs_plan() const;
    PhysicalFieldRequestPlan make_psi_det_physical_request() const;
    PhysicalFieldRequestPlan make_j_det_physical_request() const;

    void clear(Complex* out) const;
    mutable FCalculatorWorkspace workspace_;
    
    // for nonlinear term
    mutable const State* psi_nonlinear_state_ = nullptr;
    mutable double psi_nonlinear_time_ = 0.0;
    mutable bool psi_nonlinear_rhs_ready_ = false;
    mutable std::vector<Complex> psi_nonlinear_rhs_;
    mutable std::vector<double> psi_point_;

    static int temporary_field_capacity(const Params& params);

    void ensure_psi_nonlinear_rhs(const State& current, double time) const;


public:
    FCalculator(
        const Params& params,
        const Domain2D& domain,
        const SpectralMask2D& spectral_mask,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient
    );

    void rho_det(const State& current, Complex* out, double t) const;
    void psi_det(int order_parameter, const State& current, Complex* out, double t) const;
    void j_det(const State& current, Complex* out_jx, Complex* out_jy, double t) const;

    void psi_sto(int order_parameter, const State& current, Complex* out) const;
    void j_sto(const State& current, Complex* out_jx, Complex* out_jy) const;
};

#endif
