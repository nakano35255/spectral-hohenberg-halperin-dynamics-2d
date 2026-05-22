#ifndef SFI_TIME_INTEGRATOR_EULER_INCOMPRESSIBLE_H
#define SFI_TIME_INTEGRATOR_EULER_INCOMPRESSIBLE_H

#include "time_integrator.h"

class EulerIncompressible : public TimeIntegrator {
private:
    State deterministic_rhs_;
    State stochastic_rhs_;
    State u_old_;
    std::vector<Complex> rho_total_;
    std::vector<Complex> jx_increment_;
    std::vector<Complex> jy_increment_;
    bool initialized_ = false;

public:
    EulerIncompressible(
        const Domain2D& domain,
        const Params& params
    )
        : TimeIntegrator(domain, params),
          deterministic_rhs_(domain, params),
          stochastic_rhs_(domain, params),
          u_old_(domain, params),
          jx_increment_(domain.spectral_size(), Complex(0.0, 0.0)),
          jy_increment_(domain.spectral_size(), Complex(0.0, 0.0)) {}

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

        rhs.momentum_det(u_old_, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t);
        if (rhs.momentum_sto) {
            rhs.momentum_sto(u_old_, stochastic_rhs_.jx_hat_data(), stochastic_rhs_.jy_hat_data());
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

        const std::size_t jx_offset = static_cast<std::size_t>(num_components_) * local_spectral_size_;
        const std::size_t jy_offset = static_cast<std::size_t>(num_components_ + 1) * local_spectral_size_;

        for (std::size_t i = 0; i < local_spectral_size_; ++i) {
            jx_increment_[i] = dt_ * det[jx_offset + i] + sqrt_dt_ * sto[jx_offset + i];
            jy_increment_[i] = dt_ * det[jy_offset + i] + sqrt_dt_ * sto[jy_offset + i];
        }

        const Box2D& box = domain_.spectral_box();
        for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
            const double ky = domain_.ky(gy);
            for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                const double kx = domain_.kx(gx);
                const double k2 = kx * kx + ky * ky;
                if (k2 == 0.0) {
                    continue;
                }
                const std::size_t index = domain_.spectral_local_index(gx, gy);
                const Complex transverse = -ky * jx_increment_[index] + kx * jy_increment_[index];
                jx_increment_[index] = -ky * transverse / k2;
                jy_increment_[index] =  kx * transverse / k2;
            }
        }

        for (std::size_t i = 0; i < local_spectral_size_; ++i) {
            next[jx_offset + i] = current[jx_offset + i] + jx_increment_[i];
            next[jy_offset + i] = current[jy_offset + i] + jy_increment_[i];
        }

        Complex* rho_last = u.rho_hat_data(num_components_ - 1);
        std::copy(rho_total_.begin(), rho_total_.end(), rho_last);
        for (int component = 0; component < num_components_ - 1; ++component) {
            const Complex* rho = u.rho_hat_data(component);
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                rho_last[i] -= rho[i];
            }
        }

        enforce_real_symmetry(u);
    }
};

#endif
