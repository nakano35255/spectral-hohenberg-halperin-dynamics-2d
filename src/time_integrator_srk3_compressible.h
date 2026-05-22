#ifndef SFI_TIME_INTEGRATOR_SRK3_COMPRESSIBLE_H
#define SFI_TIME_INTEGRATOR_SRK3_COMPRESSIBLE_H

#include "time_integrator.h"

class SRK3Compressible : public TimeIntegrator {
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

        for (int component = 0; component < num_components_; ++component) {
            rhs.density_det(component, current, deterministic_rhs_.rho_hat_data(component), t);
        }
        rhs.momentum_det(current, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t);

        Complex* next_data = next.data();
        const Complex* current_data = current.data();
        const Complex* base_data = base.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto_a = stochastic_rhs_a_.data();
        const Complex* sto_b = stochastic_rhs_b_.data();

        for (int field = 0; field < num_fields_; ++field) {
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

        enforce_real_symmetry(next);
    }

public:
    SRK3Compressible(
        const Domain2D& domain,
        const Params& params
    )
        : TimeIntegrator(domain, params),
          stochastic_rhs_a_(domain, params),
          stochastic_rhs_b_(domain, params),
          deterministic_rhs_(domain, params),
          u_old_(domain, params),
          u_stage1_(domain, params),
          u_stage2_(domain, params) {}

    void step(State& u, double t, const RHSOperators& rhs) override {
        copy_state(u_old_, u);
        clear_state(stochastic_rhs_a_);
        clear_state(stochastic_rhs_b_);

        if (rhs.density_sto) {
            for (int component = 0; component < num_components_; ++component) {
                rhs.density_sto(component, u_old_, stochastic_rhs_a_.rho_hat_data(component));
                rhs.density_sto(component, u_old_, stochastic_rhs_b_.rho_hat_data(component));
            }
        }

        if (rhs.momentum_sto) {
            rhs.momentum_sto(u_old_, stochastic_rhs_a_.jx_hat_data(), stochastic_rhs_a_.jy_hat_data());
            rhs.momentum_sto(u_old_, stochastic_rhs_b_.jx_hat_data(), stochastic_rhs_b_.jy_hat_data());
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
