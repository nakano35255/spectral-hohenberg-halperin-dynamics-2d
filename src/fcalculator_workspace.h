#ifndef SFI_FCALCULATOR_WORKSPACE_H
#define SFI_FCALCULATOR_WORKSPACE_H

#include "domain.h"
#include "fourier_transform.h"
#include "simulationinfo.h"
#include "state.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <string>

class FCalculatorWorkspace {
private:
     const Domain2D& domain_;
     FourierTransform2D fourier_;

     int num_order_parameters_ = 0;
     int temporary_field_capacity_ = 0;

     std::size_t local_physical_size_ = 0;
     std::size_t local_spectral_size_ = 0;

     const State* cached_state_ = nullptr;
     double cached_time_ = 0.0;

     bool physical_psi_ready_ = false;

     std::vector<double> physical_psi_;
     std::vector<double> temporary_physical_;
     std::vector<Complex> temporary_spectral_;

     void reset_cache(const State& state, double time) {
          cached_state_ = &state;
          cached_time_ = time;
          physical_psi_ready_ = false;
     }

     void check_temporary_field_count(int nfields, const char* name) const {
          if (nfields < 0 || nfields > temporary_field_capacity_) {
               throw std::runtime_error(std::string("FCalculatorWorkspace ") + name + " temporary field request exceeds capacity.");
          }
     }

public:
     FCalculatorWorkspace(const Domain2D& domain, const Params& params, int temporary_field_capacity)
          : domain_(domain),
            fourier_(domain),
            num_order_parameters_(params.physics.num_order_parameters),
            temporary_field_capacity_(temporary_field_capacity),
            local_physical_size_(domain.physical_size()),
            local_spectral_size_(domain.spectral_size())
     {
          if (num_order_parameters_ < 0) {
               throw std::runtime_error("FCalculatorWorkspace requires a nonnegative number of order parameters.");
          }
          if (temporary_field_capacity_ < 0) {
               throw std::runtime_error("FCalculatorWorkspace requires nonnegative temporary field capacity.");
          }

          physical_psi_.assign(static_cast<std::size_t>(num_order_parameters_) * local_physical_size_, 0.0);

          temporary_physical_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_physical_size_, 0.0);

          temporary_spectral_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_spectral_size_, Complex(0.0, 0.0));
     }

     void ensure_cache_for(const State& state, double time) {
          if (cached_state_ == &state && cached_time_ == time) {
               return;
          }

          reset_cache(state, time);
     }

     const double* physical_psi(const State& state, double time) {
          ensure_cache_for(state, time);

          if (physical_psi_ready_) {
               return physical_psi_.data();
          }

          std::fill(physical_psi_.begin(), physical_psi_.end(), 0.0);

          if (num_order_parameters_ > 0) {
               fourier_.backward_many(num_order_parameters_, state.psi_hat_data(0), physical_psi_.data());
          }

          physical_psi_ready_ = true;
          return physical_psi_.data();
     }

     double* temporary_physical_fields(int nfields) {
          check_temporary_field_count(nfields, "physical");

          const std::size_t used_size = static_cast<std::size_t>(nfields) * local_physical_size_;

          std::fill(temporary_physical_.begin(), temporary_physical_.begin() + used_size, 0.0);

          return temporary_physical_.data();
     }

     Complex* temporary_spectral_fields(int nfields) {
          check_temporary_field_count(nfields, "spectral");

          const std::size_t used_size = static_cast<std::size_t>(nfields) * local_spectral_size_;

          std::fill(temporary_spectral_.begin(), temporary_spectral_.begin() + used_size, Complex(0.0, 0.0));

          return temporary_spectral_.data();
     }

     void forward_many(int nfields, const double* physical, Complex* spectral) {
          fourier_.forward_many(nfields, physical, spectral);
     }

     void backward_many(int nfields, const Complex* spectral, double* physical) {
          fourier_.backward_many(nfields, spectral, physical);
     }

     std::size_t local_physical_size() const {
          return local_physical_size_;
     }

     std::size_t local_spectral_size() const {
          return local_spectral_size_;
     }

     int temporary_field_capacity() const {
          return temporary_field_capacity_;
     }
};

#endif
