#include "initial_condition.h"

#include <cmath>
#include <mpi.h>
#include <stdexcept>

void validate_initial_momentum_is_transverse(const State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) {
     constexpr double ABS_TOL = 1.0e-10;
     constexpr double REL_TOL = 1.0e-10;

     const Complex* jx = state.jx_hat_data();
     const Complex* jy = state.jy_hat_data();

     int local_invalid = 0;

     for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
          if (mode.k2 == 0.0) {
               continue;
          }

          const std::size_t i = mode.index;
          const Complex k_dot_j = mode.kx * jx[i] + mode.ky * jy[i];

          const double violation = std::abs(k_dot_j);
          const double momentum_scale = std::sqrt(mode.k2) * std::sqrt(std::norm(jx[i]) + std::norm(jy[i]));
          const double tolerance = ABS_TOL + REL_TOL * momentum_scale;

          if (violation > tolerance) {
               local_invalid = 1;
               break;
          }
     }

     int global_invalid = 0;
     MPI_Allreduce(&local_invalid, &global_invalid, 1, MPI_INT, MPI_MAX, domain.comm());

     if (global_invalid) {
          throw std::runtime_error("incompressible initial momentum must be transverse: k dot j must be zero for every nonzero active spectral mode.");
     }
}