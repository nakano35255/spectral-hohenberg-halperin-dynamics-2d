#ifndef SHHD_FCALCULATOR_WORKSPACE_H
#define SHHD_FCALCULATOR_WORKSPACE_H

#include "domain.h"
#include "fcalculator_dynamics_mode.h"
#include "fcalculator_field_requests_and_layout.h"
#include "fourier_transform.h"
#include "simulationinfo.h"
#include "state.h"

#include <cstddef>
#include <vector>

class FCalculatorWorkspace {
private:
     const Domain2D& domain_;
     FourierTransform2D fourier_;

     DynamicsMode dynamics_mode_;
     PhysicalFieldLayout field_layout_;

     int num_order_parameters_ = 0;
     int temporary_field_capacity_ = 0;

     std::size_t local_physical_size_ = 0;
     std::size_t local_spectral_size_ = 0;

     const State* cached_state_ = nullptr;
     double cached_time_ = 0.0;
     bool physical_fields_ready_ = false;

     bool reference_density_ready_ = false;
     double reference_density_ = 0.0;

     std::vector<Complex> spectral_pack_;
     std::vector<double> physical_pack_;

     std::vector<double> temporary_physical_;
     std::vector<Complex> temporary_spectral_;

public:
     FCalculatorWorkspace(
          const Domain2D& domain,
          const Params& params,
          DynamicsMode dynamics_mode,
          const PhysicalRHSFieldRequests& physical_field_requests,
          int temporary_field_capacity
     );

     void ensure_physical_fields(const State& state, double time);
     double reference_density(const State& state);

     const double* physical_psi_data() const;
     const double* physical_rho_data() const;
     const double* physical_jx_data() const;
     const double* physical_jy_data() const;
     const double* physical_vx_data() const;
     const double* physical_vy_data() const;

     double* temporary_physical_fields(int nfields);
     Complex* temporary_spectral_fields(int nfields);

     void forward_many(int nfields, const double* physical, Complex* spectral);
     void backward_many(int nfields, const Complex* spectral, double* physical);

     std::size_t local_physical_size() const;
     std::size_t local_spectral_size() const;
     int temporary_field_capacity() const;
};

#endif