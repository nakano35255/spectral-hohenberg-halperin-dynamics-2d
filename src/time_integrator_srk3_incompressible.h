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
    std::vector<Complex> rho_total_;
    std::vector<Complex> jx_increment_;
    std::vector<Complex> jy_increment_;
    bool initialized_ = false;

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

        for (int component = 0; component < num_components_ - 1; ++component) {
            rhs.density_det(component, current, deterministic_rhs_.rho_hat_data(component), t);
        }
        rhs.momentum_det(current, deterministic_rhs_.jx_hat_data(), deterministic_rhs_.jy_hat_data(), t);

        Complex* next_data = next.data();
        const Complex* current_data = current.data();
        const Complex* base_data = base.data();
        const Complex* det = deterministic_rhs_.data();
        const Complex* sto_a = stochastic_rhs_a_.data();
        const Complex* sto_b = stochastic_rhs_b_.data();

        for (int field = 0; field < num_components_ - 1; ++field) {
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

        const std::size_t jx_offset = static_cast<std::size_t>(num_components_) * local_spectral_size_;
        const std::size_t jy_offset = static_cast<std::size_t>(num_components_ + 1) * local_spectral_size_;

        for (std::size_t i = 0; i < local_spectral_size_; ++i) {
            const Complex sto_jx = sto_a[jx_offset + i] + beta * sto_b[jx_offset + i];
            const Complex sto_jy = sto_a[jy_offset + i] + beta * sto_b[jy_offset + i];
            jx_increment_[i] = dt_ * det[jx_offset + i] + sqrt_dt_ * sto_jx;
            jy_increment_[i] = dt_ * det[jy_offset + i] + sqrt_dt_ * sto_jy;
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
            next_data[jx_offset + i] =
                weight_base * base_data[jx_offset + i]
              + weight_current * (current_data[jx_offset + i] + jx_increment_[i]);
            next_data[jy_offset + i] =
                weight_base * base_data[jy_offset + i]
              + weight_current * (current_data[jy_offset + i] + jy_increment_[i]);
        }

        Complex* rho_last = next.rho_hat_data(num_components_ - 1);
        std::copy(rho_total_.begin(), rho_total_.end(), rho_last);
        for (int component = 0; component < num_components_ - 1; ++component) {
            const Complex* rho = next.rho_hat_data(component);
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                rho_last[i] -= rho[i];
            }
        }

        enforce_real_symmetry(next);
    }

public:
    SRK3Incompressible(
        const Domain2D& domain,
        const Params& params
    )
        : TimeIntegrator(domain, params),
          stochastic_rhs_a_(domain, params),
          stochastic_rhs_b_(domain, params),
          deterministic_rhs_(domain, params),
          u_old_(domain, params),
          u_stage1_(domain, params),
          u_stage2_(domain, params),
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
        clear_state(stochastic_rhs_a_);
        clear_state(stochastic_rhs_b_);

        if (rhs.density_sto) {
            for (int component = 0; component < num_components_ - 1; ++component) {
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
