#ifndef SHHD_TIME_INTEGRATOR_EULER_INCOMPRESSIBLE_H
#define SHHD_TIME_INTEGRATOR_EULER_INCOMPRESSIBLE_H

#include "time_integrator.h"

class EulerIncompressible : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;
    std::vector<Complex> jx_increment_;
    std::vector<Complex> jy_increment_;

public:
    EulerIncompressible(
        const Domain2D& domain,
        const Params& params,
        const SpectralMask2D& spectral_mask
    )
        : TimeIntegrator(domain, params, spectral_mask),
          deterministic_rhs_(params, domain),
          stochastic_rhs_(params, domain),
          u_old_(params, domain),
          jx_increment_(domain.spectral_size(), Complex(0.0, 0.0)),
          jy_increment_(domain.spectral_size(), Complex(0.0, 0.0)) {}

    void step(State& u, double t, const RHSOperators& rhs) override {
        copy_state(u_old_, u);

        // optional flux output
        begin_flux_step();
        FluxBuffer* flux_det = nullptr;
        FluxBuffer* flux_sto = nullptr;
        if (flux_buffer_.requested()) {
            deterministic_flux_.begin_step();
            stochastic_flux_a_.begin_step();

            flux_det = &deterministic_flux_;
            flux_sto = &stochastic_flux_a_;
        }

        clear_state(deterministic_rhs_);
        clear_state(stochastic_rhs_);


        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            rhs.psi_det(order_parameter, u_old_, deterministic_rhs_.psi_hat_data(order_parameter), t, flux_det);
            if (rhs.psi_sto) {
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_.psi_hat_data(order_parameter), flux_sto);
            }
        }

        rhs.j_det(u_old_, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t, flux_det);
        if (rhs.j_sto) {
            rhs.j_sto(u_old_, stochastic_rhs_.jx_hat_data(), stochastic_rhs_.jy_hat_data(), flux_sto);
        }

        Complex* next = u.data();
        const Complex* current = u_old_.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto = stochastic_rhs_.data();

        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            const std::size_t offset = static_cast<std::size_t>(order_parameter + 1) * local_spectral_size_;

            for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
                const std::size_t index = offset + mode.index;
                next[index] = current[index] + dt_ * det[index] + sqrt_dt_ * sto[index];
            }
        }

        const std::size_t jx_offset = static_cast<std::size_t>(num_order_parameters_ + 1) * local_spectral_size_;
        const std::size_t jy_offset = static_cast<std::size_t>(num_order_parameters_ + 2) * local_spectral_size_;

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t index = mode.index;

            jx_increment_[index] = dt_ * det[jx_offset + index] + sqrt_dt_ * sto[jx_offset + index];
            jy_increment_[index] = dt_ * det[jy_offset + index] + sqrt_dt_ * sto[jy_offset + index];

            if (mode.k2 != 0.0) {
                const Complex transverse = -mode.ky * jx_increment_[index] + mode.kx * jy_increment_[index];

                jx_increment_[index] = -mode.ky * transverse / mode.k2;
                jy_increment_[index] =  mode.kx * transverse / mode.k2;
            }
        }

        for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
            const std::size_t index = mode.index;

            next[jx_offset + index] = current[jx_offset + index] + jx_increment_[index];
            next[jy_offset + index] = current[jy_offset + index] + jy_increment_[index];
        }

        enforce_real_symmetry(u);

        // optional flux output
        if (flux_buffer_.requested()) {
            flux_buffer_.accumulate_stage_flux(1.0, deterministic_flux_);
            flux_buffer_.accumulate_stage_flux(1.0 / sqrt_dt_, stochastic_flux_a_);
            flux_buffer_.add_incompressible_projection_pressure_flux(spectral_mask_);
        }
    }
};

#endif
