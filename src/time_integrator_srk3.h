#ifndef SFI_TIME_INTEGRATOR_SRK3_H
#define SFI_TIME_INTEGRATOR_SRK3_H

#include "time_integrator.h"

// Garcia's SRK3 coefficients.
namespace Garcia {
    constexpr double SQRT2 = 1.41421356237309504880;
    constexpr double SQRT3 = 1.73205080756887729352;

    constexpr double beta1 = ( 2.0 * SQRT2 + SQRT3)       / 5.0;
    constexpr double beta2 = (-4.0 * SQRT2 + 3.0 * SQRT3) / 5.0;
    constexpr double beta3 = ( SQRT2 - 2.0 * SQRT3)       / 10.0;
}

class SRK3Integrator : public TimeIntegrator {
private:
    State stochastic_rhs_a_;
    State stochastic_rhs_b_;
    State deterministic_rhs_;
    State u_old_;
    State u_stage1_;
    State u_stage2_;

    void calculate_stage(
        State& next,
        const State& current,
        const State& base,
        double t,
        const DetRHSFunc& calc_det,
        double beta,
        double weight_base,
        double weight_current
    ) {
        clear_state(deterministic_rhs_);
        calc_det(current, deterministic_rhs_, t);

        Complex* next_data = next.data();
        const Complex* current_data = current.data();
        const Complex* base_data = base.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto_a = stochastic_rhs_a_.data();
        const Complex* sto_b = stochastic_rhs_b_.data();

        for (int field = 0; field < num_fields_; ++field) {
            if (!mask_.active_field(field)) {
                continue;
            }

            const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                const std::size_t index = offset + i;
                const Complex stochastic = sto_a[index] + beta * sto_b[index];
                next_data[index] = weight_base * base_data[index]
                                 + weight_current * (
                                       current_data[index] + dt_ * det[index] + sqrt_dt_ * stochastic
                                   );
            }
        }

        enforce_after_stage(next);
    }

public:
    SRK3Integrator(
        const Domain2D& domain,
        const Params& params,
        const FluidConstraint& constraint,
        const TimeEvolutionMask& mask
    )
        : TimeIntegrator(domain, params, constraint, mask),
          stochastic_rhs_a_(domain, params),
          stochastic_rhs_b_(domain, params),
          deterministic_rhs_(domain, params),
          u_old_(domain, params),
          u_stage1_(domain, params),
          u_stage2_(domain, params) {}

    void step(State& u, double t, const DetRHSFunc& calc_det, const StoRHSFunc& calc_sto) override {
        copy_state(u_old_, u);
        clear_state(stochastic_rhs_a_);
        clear_state(stochastic_rhs_b_);

        if (calc_sto) {
            calc_sto(stochastic_rhs_a_);
            calc_sto(stochastic_rhs_b_);
        }

        calculate_stage(
            u_stage1_, u_old_, u_old_, t,
            calc_det, Garcia::beta1, 0.0, 1.0
        );
        calculate_stage(
            u_stage2_, u_stage1_, u_old_, t + dt_,
            calc_det, Garcia::beta2, 3.0 / 4.0, 1.0 / 4.0
        );
        calculate_stage(
            u, u_stage2_, u_old_, t + 0.5 * dt_,
            calc_det, Garcia::beta3, 1.0 / 3.0, 2.0 / 3.0
        );
    }
};

#endif
