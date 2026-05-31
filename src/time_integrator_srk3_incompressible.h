#ifndef SHHD_TIME_INTEGRATOR_SRK3_INCOMPRESSIBLE_H
#define SHHD_TIME_INTEGRATOR_SRK3_INCOMPRESSIBLE_H

#include "time_integrator.h"

class SRK3Incompressible : public TimeIntegrator {
private:
    static constexpr double SQRT2 = 1.41421356237309504880;
    static constexpr double SQRT3 = 1.73205080756887729352;
    static constexpr double BETA1 = ( 2.0 * SQRT2 + SQRT3)       / 5.0;
    static constexpr double BETA2 = (-4.0 * SQRT2 + 3.0 * SQRT3) / 5.0;
    static constexpr double BETA3 = ( SQRT2 - 2.0 * SQRT3)       / 10.0;
    static constexpr double STAGE1_BASE_WEIGHT = 0.0;
    static constexpr double STAGE1_CURRENT_WEIGHT = 1.0;
    static constexpr double STAGE1_FLUX_WEIGHT = 1.0 / 6.0;
    static constexpr double STAGE2_BASE_WEIGHT = 3.0 / 4.0;
    static constexpr double STAGE2_CURRENT_WEIGHT = 1.0 / 4.0;
    static constexpr double STAGE2_FLUX_WEIGHT = 1.0 / 6.0;
    static constexpr double STAGE3_BASE_WEIGHT = 1.0 / 3.0;
    static constexpr double STAGE3_CURRENT_WEIGHT = 2.0 / 3.0;
    static constexpr double STAGE3_FLUX_WEIGHT = 2.0 / 3.0;

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
        double weight_current,
        double flux_weight
    ) {
        clear_state(deterministic_rhs_);
        copy_state(next, current);

        // optional flux output
        FluxBuffer* flux_det = nullptr;
        if (flux_buffer_.requested()) {
            deterministic_flux_.begin_step();
            flux_det = &deterministic_flux_;
        }

        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            rhs.psi_det(order_parameter, current, deterministic_rhs_.psi_hat_data(order_parameter), t, flux_det);
        }
        rhs.j_det(current, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t, flux_det);

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

        // optional flux output
        if (flux_buffer_.requested()) {
            flux_buffer_.accumulate_stage_flux(flux_weight, deterministic_flux_);
        }
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

        // optional flux output
        begin_flux_step();

        clear_state(stochastic_rhs_a_);
        clear_state(stochastic_rhs_b_);

        // optional flux output
        FluxBuffer* flux_sto_a = nullptr;
        if (flux_buffer_.requested()) {
            stochastic_flux_a_.begin_step();

            flux_sto_a = &stochastic_flux_a_;
        }

        if (rhs.psi_sto) {
            for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_a_.psi_hat_data(order_parameter), flux_sto_a);
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_b_.psi_hat_data(order_parameter), nullptr);
            }
        }

        if (rhs.j_sto) {
            rhs.j_sto(u_old_, stochastic_rhs_a_.jx_hat_data(), stochastic_rhs_a_.jy_hat_data(), flux_sto_a);
            rhs.j_sto(u_old_, stochastic_rhs_b_.jx_hat_data(), stochastic_rhs_b_.jy_hat_data(), nullptr);
        }

        calculate_stage(
            u_stage1_, u_old_, u_old_, t,
            rhs, BETA1, STAGE1_BASE_WEIGHT, STAGE1_CURRENT_WEIGHT, STAGE1_FLUX_WEIGHT
        );

        calculate_stage(
            u_stage2_, u_stage1_, u_old_, t + dt_,
            rhs, BETA2, STAGE2_BASE_WEIGHT, STAGE2_CURRENT_WEIGHT, STAGE2_FLUX_WEIGHT
        );

        calculate_stage(
            u, u_stage2_, u_old_, t + 0.5 * dt_,
            rhs, BETA3, STAGE3_BASE_WEIGHT, STAGE3_CURRENT_WEIGHT, STAGE3_FLUX_WEIGHT
        );

        // optional flux output
        if (flux_buffer_.requested()) {
            flux_buffer_.accumulate_stage_flux(1.0 / sqrt_dt_, stochastic_flux_a_);
            flux_buffer_.add_incompressible_projection_pressure_flux(spectral_mask_);
        }
    }
};

#endif
