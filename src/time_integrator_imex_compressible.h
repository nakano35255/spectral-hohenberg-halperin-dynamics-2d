#ifndef SHHD_TIME_INTEGRATOR_IMEX_COMPRESSIBLE_H
#define SHHD_TIME_INTEGRATOR_IMEX_COMPRESSIBLE_H

#include "time_integrator.h"

class IMEXCompressible : public TimeIntegrator {
private:
     State stochastic_rhs_a_;
     State stochastic_rhs_b_;
     State u_old_;
     State u_half_;
     State nonlinear_rhs_old_;
     State nonlinear_rhs_half_;
     State linear_rhs_old_;
     State implicit_rhs_;

     // optional flux output
     State linear_rhs_next_;
     FluxBuffer stochastic_flux_b_;

     void calculate_nonlinear_rhs(State& out, const State& current, double t, const RHSOperators& rhs, FluxBuffer* flux) {
          clear_state(out);

          rhs.rho_nonlin_det(current, out.rho_hat_data(), t, flux);

          for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
               rhs.psi_nonlin_det(order_parameter, current, out.psi_hat_data(order_parameter), t, flux);
          }

          rhs.j_nonlin_det(current, out.jx_hat_data(), out.jy_hat_data(), t, flux);
     }

     void calculate_linear_rhs(State& out, const State& current, double t, const RHSOperators& rhs, FluxBuffer* flux) {
          clear_state(out);

          rhs.rho_lin_det(current, out.rho_hat_data(), t, flux);

          for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
               rhs.psi_lin_det(order_parameter, current, out.psi_hat_data(order_parameter), t, flux);
          }

          rhs.j_lin_det(current, out.jx_hat_data(), out.jy_hat_data(), t, flux);
     }

     void calculate_stochastic_rhs(const State& current, const RHSOperators& rhs) {
          clear_state(stochastic_rhs_a_);
          clear_state(stochastic_rhs_b_);

          // optional flux output
          FluxBuffer* flux_sto_a = nullptr;
          FluxBuffer* flux_sto_b = nullptr;
          if (flux_buffer_.requested()) {
               stochastic_flux_a_.begin_step();
               stochastic_flux_b_.set_request(flux_buffer_.request());
               stochastic_flux_b_.begin_step();

               flux_sto_a = &stochastic_flux_a_;
               flux_sto_b = &stochastic_flux_b_;
          }

          if (rhs.psi_sto) {
               for (int order_parameter = 0; order_parameter < num_order_parameters_; ++order_parameter) {
                    rhs.psi_sto(order_parameter, current, stochastic_rhs_a_.psi_hat_data(order_parameter), flux_sto_a);
                    rhs.psi_sto(order_parameter, current, stochastic_rhs_b_.psi_hat_data(order_parameter), flux_sto_b);
               }
          }

          if (rhs.j_sto) {
               rhs.j_sto(current, stochastic_rhs_a_.jx_hat_data(), stochastic_rhs_a_.jy_hat_data(), flux_sto_a);
               rhs.j_sto(current, stochastic_rhs_b_.jx_hat_data(), stochastic_rhs_b_.jy_hat_data(), flux_sto_b);
          }
     }

     void predictor_step(State& u_half, double t, const RHSOperators& rhs) {
          const double alpha = 0.5 * dt_;
          const double sqrt_alpha = std::sqrt(alpha);

          calculate_nonlinear_rhs(nonlinear_rhs_old_, u_old_, t, rhs, nullptr);

          clear_state(implicit_rhs_);

          Complex* out = implicit_rhs_.data();
          const Complex* old_data = u_old_.data();
          const Complex* nonlinear = nonlinear_rhs_old_.data();
          const Complex* noise = stochastic_rhs_a_.data();

          for (int field = 0; field < num_fields_; ++field) {
               const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;

               for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
                    const std::size_t index = offset + mode.index;
                    out[index] = old_data[index]
                              + alpha * nonlinear[index]
                              + sqrt_alpha * noise[index];
               }
          }

          enforce_real_symmetry(implicit_rhs_);

          rhs.apply_implicit_inverse(u_half, implicit_rhs_, alpha);
          enforce_real_symmetry(u_half);
     }

     void corrector_step(State& u_next, double t, const RHSOperators& rhs) {
          const double alpha = 0.5 * dt_;
          const double sqrt_alpha = std::sqrt(alpha);

          // optional flux output
          FluxBuffer* flux = nullptr;
          if (flux_buffer_.requested()) {
               deterministic_flux_.begin_step();
               flux = &deterministic_flux_;
          }

          calculate_linear_rhs(linear_rhs_old_, u_old_, t, rhs, flux);
          // optional flux output
          if (flux_buffer_.requested()) {
               flux_buffer_.accumulate_stage_flux(0.5, deterministic_flux_);
               deterministic_flux_.begin_step();
          }

          calculate_nonlinear_rhs(nonlinear_rhs_half_, u_half_, t + alpha, rhs, flux);
          // optional flux output
          if (flux_buffer_.requested()) {
               flux_buffer_.accumulate_stage_flux(1.0, deterministic_flux_);
          }

          clear_state(implicit_rhs_);

          Complex* out = implicit_rhs_.data();
          const Complex* old_data = u_old_.data();
          const Complex* linear_old = linear_rhs_old_.data();
          const Complex* nonlinear_half = nonlinear_rhs_half_.data();
          const Complex* noise_a = stochastic_rhs_a_.data();
          const Complex* noise_b = stochastic_rhs_b_.data();

          for (int field = 0; field < num_fields_; ++field) {
               const std::size_t offset = static_cast<std::size_t>(field) * local_spectral_size_;

               for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
                    const std::size_t index = offset + mode.index;
                    out[index] = old_data[index]
                              + alpha * linear_old[index]
                              + dt_ * nonlinear_half[index]
                              + sqrt_alpha * (noise_a[index] + noise_b[index]);
               }
          }

          enforce_real_symmetry(implicit_rhs_);

          rhs.apply_implicit_inverse(u_next, implicit_rhs_, alpha);
          enforce_real_symmetry(u_next);

          // optional flux output
          if (flux_buffer_.requested()) {
               deterministic_flux_.begin_step();
               calculate_linear_rhs(linear_rhs_next_, u_next, t + dt_, rhs, &deterministic_flux_);
               flux_buffer_.accumulate_stage_flux(0.5, deterministic_flux_);
          }
     }

public:
     IMEXCompressible(
          const Domain2D& domain,
          const Params& params,
          const SpectralMask2D& spectral_mask
     )
          : TimeIntegrator(domain, params, spectral_mask),
          stochastic_rhs_a_(params, domain),
          stochastic_rhs_b_(params, domain),
          u_old_(params, domain),
          u_half_(params, domain),
          nonlinear_rhs_old_(params, domain),
          nonlinear_rhs_half_(params, domain),
          linear_rhs_old_(params, domain),
          implicit_rhs_(params, domain),
          linear_rhs_next_(params, domain),
          stochastic_flux_b_(params, domain) {}

     void step(State& u, double t, const RHSOperators& rhs) override {          
          // optional flux output
          begin_flux_step();

          copy_state(u_old_, u);

          calculate_stochastic_rhs(u_old_, rhs);

          predictor_step(u_half_, t, rhs);
          corrector_step(u, t, rhs);

          // optional flux output
          if (flux_buffer_.requested()) {
               flux_buffer_.accumulate_stage_flux(std::sqrt(0.5) / sqrt_dt_, stochastic_flux_a_);
               flux_buffer_.accumulate_stage_flux(std::sqrt(0.5) / sqrt_dt_, stochastic_flux_b_);
          }
     }
};

#endif