#ifndef SFI_FCALCULATOR_WORKSPACE_H
#define SFI_FCALCULATOR_WORKSPACE_H

#include "domain.h"
#include "fourier_transform.h"
#include "simulationinfo.h"
#include "state.h"

#include <cstddef>
#include <vector>

struct PhysicalFieldRequestPlan {
     bool need_psi = false;
     bool need_rho = false;
     bool need_j = false;
     bool need_velocity = false;

     bool any() const {
          return need_psi || need_rho || need_j || need_velocity;
     }
};

struct PhysicalRHSPlan {
     PhysicalFieldRequestPlan psi_det;
     PhysicalFieldRequestPlan j_det;
};

struct PhysicalFieldOffsets {
     int psi = -1;
     int rho = -1;
     int jx = -1;
     int jy = -1;
     int vx = -1;
     int vy = -1;
};

class FCalculatorWorkspace {
private:
     const Domain2D& domain_;
     const Params& params_;
     FourierTransform2D fourier_;

     PhysicalRHSPlan physical_rhs_plan_;
     PhysicalFieldRequestPlan request_;
     PhysicalFieldOffsets offsets_;

     int num_order_parameters_ = 0;
     int temporary_field_capacity_ = 0;

     std::size_t local_physical_size_ = 0;
     std::size_t local_spectral_size_ = 0;

     int ntransform_ = 0;
     int nphysical_ = 0;

     bool is_quiescent_ = false;
     bool is_incompressible_ = false;
     bool is_compressible_ = false;

     const State* cached_state_ = nullptr;
     double cached_time_ = 0.0;
     bool physical_fields_ready_ = false;

     bool reference_density_needed_ = false;
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
          const PhysicalRHSPlan& physical_rhs_plan,
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
