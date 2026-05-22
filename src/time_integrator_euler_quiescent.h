#ifndef SFI_TIME_INTEGRATOR_EULER_QUIESCENT_H
#define SFI_TIME_INTEGRATOR_EULER_QUIESCENT_H

#include "time_integrator.h"

class EulerQuiescent : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;
    std::vector<Complex> rho_total_;
    bool initialized_ = false;

public:
    EulerQuiescent(
        const Domain2D& domain,
        const Params& params
    )
        : TimeIntegrator(domain, params),
          deterministic_rhs_(domain, params),
          stochastic_rhs_(domain, params),
          u_old_(domain, params) {}

    void step(State& u, double t, const RHSOperators& rhs) override {
        if (!initialized_) {
            rho_total_.assign(local_spectral_size_, Complex(0.0, 0.0));
            for (int component = 0; component < num_components_; ++component) {
                const Complex* rho = u.rho_hat_data(component);
                for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                    rho_total_[i] += rho[i];
                }
            }
            initialized_ = true;
        }

        copy_state(u_old_, u);
        clear_state(deterministic_rhs_);
        clear_state(stochastic_rhs_);

        for (int component = 0; component < num_components_ - 1; ++component) {
            rhs.density_det(component, u_old_, deterministic_rhs_.rho_hat_data(component), t);
            if (rhs.density_sto) {
                rhs.density_sto(component, u_old_, stochastic_rhs_.rho_hat_data(component));
            }
        }

        Complex* next = u.data();
        const Complex* current = u_old_.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto = stochastic_rhs_.data();

        for (int field = 0; field < num_components_ - 1; ++field) {
            const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                const std::size_t index = offset + i;
                next[index] = current[index] + dt_ * det[index] + sqrt_dt_ * sto[index];
            }
        }

        Complex* rho_last = u.rho_hat_data(num_components_ - 1);
        std::copy(rho_total_.begin(), rho_total_.end(), rho_last);
        for (int component = 0; component < num_components_ - 1; ++component) {
            const Complex* rho = u.rho_hat_data(component);
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                rho_last[i] -= rho[i];
            }
        }

        std::fill(u.jx_hat_data(), u.jx_hat_data() + local_spectral_size_, Complex(0.0, 0.0));
        std::fill(u.jy_hat_data(), u.jy_hat_data() + local_spectral_size_, Complex(0.0, 0.0));

        enforce_real_symmetry(u);
    }
};

#endif
