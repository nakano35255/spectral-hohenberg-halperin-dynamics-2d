#include "implicit_operator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------- //
double ImplicitOperator::reference_density(const State& state) const {
     if (reference_density_ready_) {
          return reference_density_;
     }

     const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

     double local_rho0 = 0.0;

     const Box2D& box = domain_.spectral_box();
     if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
          const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
          const std::size_t index = static_cast<std::size_t>(0 - box.low[1]) * local_nkx + static_cast<std::size_t>(0 - box.low[0]);

          local_rho0 = state.rho_hat_data()[index].real() / grid_size;
     }

     MPI_Allreduce(&local_rho0, &reference_density_, 1, MPI_DOUBLE, MPI_SUM, domain_.comm());

     if (reference_density_ == 0.0) {
          throw std::runtime_error("ImplicitOperator requires nonzero reference density.");
     }

     reference_density_ready_ = true;
     return reference_density_;
}
// ---------------------------------------------------------------------- //
void ImplicitOperator::solve_order_parameter_fields(State& dst, const State& src, double alpha) const {
     const std::vector<double>& mobility = transport_coefficient_.order_parameter_mobility();

     for (int order_parameter = 0; order_parameter < src.num_order_parameters(); ++order_parameter) {
          const double mobility_op = mobility[static_cast<std::size_t>(order_parameter)];
          const double mu_k0 = free_energy_.chemical_potential_k0_coefficient(order_parameter);
          const double mu_k2 = free_energy_.chemical_potential_k2_coefficient(order_parameter);
          const double mu_k4 = free_energy_.chemical_potential_k4_coefficient(order_parameter);

          Complex* psi_dst = dst.psi_hat_data(order_parameter);
          const Complex* psi_src = src.psi_hat_data(order_parameter);

          if (mobility_op == 0.0 || (mu_k0 == 0.0 && mu_k2 == 0.0 && mu_k4 == 0.0)) {
               continue;
          }

          for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
               const std::size_t i = mode.index;
               const double k2 = mode.k2;
               const double chemical = mu_k0 + mu_k2 * k2 + mu_k4 * k2 * k2;

               const double denominator = 1.0 + alpha * mobility_op * k2 * chemical;

               psi_dst[i] = psi_src[i] / denominator;
          }
     }
}
// ---------------------------------------------------------------------- //
void ImplicitOperator::solve_incompressible_hydrodynamic_fields(State& dst, const State& src, double alpha) const {
     const double rho0 = reference_density(src);
     const double eta = transport_coefficient_.shear_viscosity();
     const double nu = eta / rho0;

     Complex* jx_dst = dst.jx_hat_data();
     Complex* jy_dst = dst.jy_hat_data();

     const Complex* jx_src = src.jx_hat_data();
     const Complex* jy_src = src.jy_hat_data();

     for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
          const std::size_t i = mode.index;

          if (mode.k2 == 0.0) {
               jx_dst[i] = jx_src[i];
               jy_dst[i] = jy_src[i];
               continue;
          }

          const Complex transverse = -mode.ky * jx_src[i] + mode.kx * jy_src[i];
          const double denominator = 1.0 + alpha * nu * mode.k2;

          jx_dst[i] = (-mode.ky * transverse / mode.k2) / denominator;
          jy_dst[i] = ( mode.kx * transverse / mode.k2) / denominator;
     }
}
// ---------------------------------------------------------------------- //
void ImplicitOperator::solve_compressible_hydrodynamic_fields(State& dst, const State& src, double alpha) const {
     const double rho0 = reference_density(src);
     const double c2 = thermodynamics_.linear_pressure_coefficient();
     const double eta = transport_coefficient_.shear_viscosity();
     const double zeta = transport_coefficient_.bulk_viscosity();

     const double nu_t = eta / rho0;
     const double nu_l = (eta + zeta) / rho0;

     Complex* rho_dst = dst.rho_hat_data();
     Complex* jx_dst = dst.jx_hat_data();
     Complex* jy_dst = dst.jy_hat_data();

     const Complex* rho_src = src.rho_hat_data();
     const Complex* jx_src = src.jx_hat_data();
     const Complex* jy_src = src.jy_hat_data();

     for (const SpectralMode2D& mode : spectral_mask_.active_modes()) {
          const std::size_t i = mode.index;
          const double k2 = mode.k2;

          if (k2 == 0.0) {
               rho_dst[i] = rho_src[i];
               jx_dst[i] = jx_src[i];
               jy_dst[i] = jy_src[i];
               continue;
          }

          const double k = std::sqrt(k2);
          const double ex = mode.kx / k;
          const double ey = mode.ky / k;

          const Complex b_rho = rho_src[i];
          const Complex b_jx = jx_src[i];
          const Complex b_jy = jy_src[i];

          const Complex b_q = ex * b_jx + ey * b_jy;
          const Complex b_s = -ey * b_jx + ex * b_jy;

          const double a = 1.0 + alpha * nu_l * k2;
          const double determinant = a + alpha * alpha * c2 * k2;
          const double transverse_denominator = 1.0 + alpha * nu_t * k2;

          const Complex imaginary(0.0, 1.0);

          const Complex rho = (a * b_rho - imaginary * alpha * k * b_q) / determinant;
          const Complex q = (b_q - imaginary * alpha * c2 * k * b_rho) / determinant;
          const Complex s = b_s / transverse_denominator;

          rho_dst[i] = rho;
          jx_dst[i] = ex * q - ey * s;
          jy_dst[i] = ey * q + ex * s;
     }
}
// ---------------------------------------------------------------------- //
void ImplicitOperator::apply_inverse(State& dst, const State& src, double alpha) const {
     const std::size_t total_size = static_cast<std::size_t>(src.num_fields()) * domain_.spectral_size();

     std::copy(src.data(), src.data() + total_size, dst.data());

     solve_order_parameter_fields(dst, src, alpha);

     if (is_compressible_mode(dynamics_mode_)) {
          solve_compressible_hydrodynamic_fields(dst, src, alpha);
     }
     else if (is_incompressible_mode(dynamics_mode_)) {
          solve_incompressible_hydrodynamic_fields(dst, src, alpha);
     }
     else if (is_quiescent_mode(dynamics_mode_)) {
          // Nothing else to solve.
     }
}