#ifndef SFI_TIME_INTEGRATOR_SRK3_INCOMPRESSIBLE_H
#define SFI_TIME_INTEGRATOR_SRK3_INCOMPRESSIBLE_H

#include "time_integrator.h"

class SRK3Incompressible : public TimeIntegrator {
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
    std::vector<Complex> jx_increment_;
    std::vector<Complex> jy_increment_;

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
        rhs.j_det(current, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t);

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

        const std::size_t jx_offset = static_cast<std::size_t>(num_order_parameters_ + 1) * local_spectral_size_;
        const std::size_t jy_offset = static_cast<std::size_t>(num_order_parameters_ + 2) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t index = mode.index;

            const Complex sto_jx = sto_a[jx_offset + index] + beta * sto_b[jx_offset + index];
            const Complex sto_jy = sto_a[jy_offset + index] + beta * sto_b[jy_offset + index];

            jx_increment_[index] = dt_ * det[jx_offset + index] + sqrt_dt_ * sto_jx;
            jy_increment_[index] = dt_ * det[jy_offset + index] + sqrt_dt_ * sto_jy;

            if (mode.k2 != 0.0) {
                const Complex transverse = -mode.ky * jx_increment_[index] + mode.kx * jy_increment_[index];

                jx_increment_[index] = -mode.ky * transverse / mode.k2;
                jy_increment_[index] =  mode.kx * transverse / mode.k2;
            }
        }

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t index = mode.index;

            next_data[jx_offset + index] =
                weight_base * base_data[jx_offset + index]
              + weight_current * (current_data[jx_offset + index] + jx_increment_[index]);
            next_data[jy_offset + index] =
                weight_base * base_data[jy_offset + index]
              + weight_current * (current_data[jy_offset + index] + jy_increment_[index]);
        }

        enforce_real_symmetry(next);
    }

public:
    SRK3Incompressible(
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
          u_stage2_(params, domain),
          jx_increment_(domain.spectral_size(), Complex(0.0, 0.0)),
          jy_increment_(domain.spectral_size(), Complex(0.0, 0.0)) {}

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

        if (rhs.j_sto) {
            rhs.j_sto(u_old_, stochastic_rhs_a_.jx_hat_data(), stochastic_rhs_a_.jy_hat_data());
            rhs.j_sto(u_old_, stochastic_rhs_b_.jx_hat_data(), stochastic_rhs_b_.jy_hat_data());
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
