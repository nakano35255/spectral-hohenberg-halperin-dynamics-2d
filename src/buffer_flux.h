#ifndef SHHD_BUFFER_FLUX_H
#define SHHD_BUFFER_FLUX_H

#include "domain.h"
#include "simulationinfo.h"
#include "spectral_mask.h"

#include <cstddef>
#include <vector>

struct FluxRequest {
     bool density = false;
     bool order_parameter = false;
     bool momentum = false;

     bool any() const;
};

class FluxBuffer {
private:
     int num_order_parameters_ = 0;
     std::size_t local_spectral_size_ = 0;

     FluxRequest request_;
     std::vector<Complex> data_;

     int density_x_slot() const;
     int density_y_slot() const;
     int order_parameter_x_slot(int order_parameter) const;
     int order_parameter_y_slot(int order_parameter) const;
     int momentum_base_slot() const;
     int momentum_xx_slot() const;
     int momentum_xy_slot() const;
     int momentum_yy_slot() const;

     int num_slots() const;
     std::size_t offset(int slot) const;

     void check_order_parameter(int order_parameter) const;
     void check_compatible(const FluxBuffer& flux) const;

     const Complex* slot_data(int slot) const;
     Complex* slot_data(int slot);

     void accumulate_slot(int slot, double weight, const Complex* values);

public:
     FluxBuffer(const Params& params, const Domain2D& domain);

     void set_request(const FluxRequest& request);
     const FluxRequest& request() const;
     bool requested() const;

     void begin_step();

     Complex* density_flux_x_hat_data();
     const Complex* density_flux_x_hat_data() const;
     Complex* density_flux_y_hat_data();
     const Complex* density_flux_y_hat_data() const;
     Complex* order_parameter_flux_x_hat_data(int order_parameter);
     const Complex* order_parameter_flux_x_hat_data(int order_parameter) const;
     Complex* order_parameter_flux_y_hat_data(int order_parameter);
     const Complex* order_parameter_flux_y_hat_data(int order_parameter) const;
     Complex* momentum_flux_xx_hat_data();
     const Complex* momentum_flux_xx_hat_data() const;
     Complex* momentum_flux_xy_hat_data();
     const Complex* momentum_flux_xy_hat_data() const;
     Complex* momentum_flux_yy_hat_data();
     const Complex* momentum_flux_yy_hat_data() const;

     void accumulate_stage_flux(double weight, const FluxBuffer& stage_flux);

     void add_incompressible_projection_pressure_flux(const SpectralMask2D& spectral_mask);
};

#endif
