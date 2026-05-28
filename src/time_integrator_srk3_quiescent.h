#ifndef SFI_TIME_INTEGRATOR_SRK3_QUIESCENT_H
#define SFI_TIME_INTEGRATOR_SRK3_QUIESCENT_H

#include "time_integrator.h"

class SRK3Quiescent : public TimeIntegrator {
private:
    static constexpr double SQRT2 = 1.41421356237309504880;
    static constexpr double SQRT3 = 1.73205080756887729352;
    static constexpr double BETA1 = ( 2.0 * SQRT2 + SQRT3)       / 5.0;
    static constexpr double BETA2 = (-4.0 * SQRT2 + 3.0 * SQRT3) / 5.0;
    static constexpr double BETA3 = ( SQRT2 - 2.0 * SQRT3)       / 10.0;

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
        const RHSOperators& rhs,
        double beta,
        double weight_base,
        double weight_current
    ) {
        clear_state(deterministic_rhs_);
        copy_state(next, current);

        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            rhs.psi_det(order_parameter, current, deterministic_rhs_.psi_hat_data(order_parameter), t);
        }

        Complex* next_data = next.data();
        const Complex* current_data = current.data();
        const Complex* base_data = base.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto_a = stochastic_rhs_a_.data();
        const Complex* sto_b = stochastic_rhs_b_.data();

        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            const std::size_t offset = static_cast<std::size_t>(order_parameter + 1) * local_spectral_size_;

            for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
                const std::size_t index = offset + mode.index;
                const Complex stochastic = sto_a[index] + beta * sto_b[index];
                next_data[index] = weight_base * base_data[index]
                                 + weight_current * (
                                       current_data[index] + dt_ * det[index] + sqrt_dt_ * stochastic
                                   );
            }
        }

        enforce_real_symmetry(next);
    }

public:
    SRK3Quiescent(
        const Domain2D& domain,
        const Params& params,
        const SpectralMask2D& spectral_mask
    )
        : TimeIntegrator(domain, params, spectral_mask),
          stochastic_rhs_a_(params, domain),
          stochastic_rhs_b_(params, domain),
          deterministic_rhs_(params, domain),
          u_old_(params, domain),
          u_stage1_(params, domain),
          u_stage2_(params, domain) {}

    void step(State& u, double t, const RHSOperators& rhs) override {
        copy_state(u_old_, u);
        clear_state(stochastic_rhs_a_);
        clear_state(stochastic_rhs_b_);

        if (rhs.psi_sto) {
            for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_a_.psi_hat_data(order_parameter));
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_b_.psi_hat_data(order_parameter));
            }
        }

        calculate_stage(
            u_stage1_, u_old_, u_old_, t,
            rhs, BETA1, 0.0, 1.0
        );
        calculate_stage(
            u_stage2_, u_stage1_, u_old_, t + dt_,
            rhs, BETA2, 3.0 / 4.0, 1.0 / 4.0
        );
        calculate_stage(
            u, u_stage2_, u_old_, t + 0.5 * dt_,
            rhs, BETA3, 1.0 / 3.0, 2.0 / 3.0
        );
    }
};

#endif
