#ifndef SHHD_TIME_INTEGRATOR_SRK3_COMPRESSIBLE_H
#define SHHD_TIME_INTEGRATOR_SRK3_COMPRESSIBLE_H

#include "time_integrator.h"

class SRK3Compressible : public TimeIntegrator {
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

        // optional flux output
        FluxBuffer* flux_det = nullptr;
        if (flux_buffer_.requested()) {
            deterministic_flux_.begin_step();
            flux_det = &deterministic_flux_;
        }

        rhs.rho_det(current, deterministic_rhs_.rho_hat_data(), t, flux_det);

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

        for (int field = 0; field < num_fields_; ++field) {
            const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;

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

        // optional flux output
        if (flux_buffer_.requested()) {
            flux_buffer_.accumulate_stage_flux(flux_weight, deterministic_flux_);
        }
    }

public:
    SRK3Compressible(
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
        }
    }
};

#endif
