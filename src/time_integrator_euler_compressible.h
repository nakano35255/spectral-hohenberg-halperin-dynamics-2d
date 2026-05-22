#ifndef SFI_TIME_INTEGRATOR_EULER_COMPRESSIBLE_H
#define SFI_TIME_INTEGRATOR_EULER_COMPRESSIBLE_H

#include "time_integrator.h"

class EulerCompressible : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;

public:
    EulerCompressible(
        const Domain2D& domain,
        const Params& params
    )
        : TimeIntegrator(domain, params),
          deterministic_rhs_(domain, params),
          stochastic_rhs_(domain, params),
          u_old_(domain, params) {}

    void step(State& u, double t, const RHSOperators& rhs) override {

        copy_state(u_old_, u);
        clear_state(deterministic_rhs_);
        clear_state(stochastic_rhs_);

        for (int component = 0; component < num_components_; ++component) {
            rhs.density_det(component, u_old_, deterministic_rhs_.rho_hat_data(component), t);
            if (rhs.density_sto) {
                rhs.density_sto(component, u_old_, stochastic_rhs_.rho_hat_data(component));
            }
        }

        rhs.momentum_det(u_old_, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t);
        if (rhs.momentum_sto) {
            rhs.momentum_sto(u_old_, stochastic_rhs_.jx_hat_data(), stochastic_rhs_.jy_hat_data());
        }

        Complex* next = u.data();
        const Complex* current = u_old_.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto = stochastic_rhs_.data();

        for (int field = 0; field < num_fields_; ++field) {
            const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                const std::size_t index = offset + i;
                next[index] = current[index] + dt_ * det[index] + sqrt_dt_ * sto[index];
            }
        }

        enforce_real_symmetry(u);
    }
};

#endif
