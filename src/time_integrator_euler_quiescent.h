#ifndef SFI_TIME_INTEGRATOR_EULER_QUIESCENT_H
#define SFI_TIME_INTEGRATOR_EULER_QUIESCENT_H

#include "time_integrator.h"

class EulerQuiescent : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;

public:
    EulerQuiescent(
        const Domain2D& domain,
        const Params& params,
        const SpectralMask2D& spectral_mask
    )
        : TimeIntegrator(domain, params, spectral_mask),
          deterministic_rhs_(domain, params),
          stochastic_rhs_(domain, params),
          u_old_(domain, params) {}

    void step(State& u, double t, const RHSOperators& rhs) override {
        copy_state(u_old_, u);
        clear_state(deterministic_rhs_);
        clear_state(stochastic_rhs_);

        for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
            rhs.psi_det(order_parameter, u_old_, deterministic_rhs_.psi_hat_data(order_parameter), t);
            if (rhs.psi_sto) {
                rhs.psi_sto(order_parameter, u_old_, stochastic_rhs_.psi_hat_data(order_parameter));
            }
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

        enforce_real_symmetry(u);
    }
};

#endif
