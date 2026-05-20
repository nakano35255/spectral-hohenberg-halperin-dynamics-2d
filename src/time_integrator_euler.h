#ifndef SFI_TIME_INTEGRATOR_EULER_H
#define SFI_TIME_INTEGRATOR_EULER_H

#include "time_integrator.h"

class EulerIntegrator : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;

public:
    EulerIntegrator(
        const Domain2D& domain,
        const Params& params,
        const FluidConstraint& constraint,
        const TimeEvolutionMask& mask
    )
        : TimeIntegrator(domain, params, constraint, mask),
          deterministic_rhs_(domain, params),
          stochastic_rhs_(domain, params),
          u_old_(domain, params) {}

    void step(State& u, double t, const DetRHSFunc& calc_det, const StoRHSFunc& calc_sto) override {
        copy_state(u_old_, u);
        clear_state(deterministic_rhs_);
        clear_state(stochastic_rhs_);

        calc_det(u_old_, deterministic_rhs_, t);
        if (calc_sto) {
            calc_sto(stochastic_rhs_);
        }

        Complex* next = u.data();
        const Complex* current = u_old_.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto = stochastic_rhs_.data();

        for (int field = 0; field < num_fields_; ++field) {
            if (!mask_.active_field(field)) {
                continue;
            }

            const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                const std::size_t index = offset + i;
                next[index] = current[index] + dt_ * det[index] + sqrt_dt_ * sto[index];
            }
        }

        enforce_after_stage(u);
    }
};

#endif
