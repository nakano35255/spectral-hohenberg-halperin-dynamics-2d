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
    RHSEvaluationPlan rhs_plan_;

    void clear(Complex* out) const;

    RHSEvaluationPlan make_rhs_evaluation_plan() const;
    PhysicalFieldPlan make_physical_field_plan() const;
    bool is_quiescent_time_evolution() const;

    mutable FCalculatorWorkspace workspace_;

    mutable const State* psi_nonlinear_state_ = nullptr;
    mutable double psi_nonlinear_time_ = 0.0;
    mutable bool psi_nonlinear_rhs_ready_ = false;
    mutable std::vector<Complex> psi_nonlinear_rhs_;
    mutable std::vector<double> psi_point_;
    mutable std::uint64_t psi_noise_call_count_ = 0;

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
};

#endif
