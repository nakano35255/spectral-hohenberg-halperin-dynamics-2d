#include "buffer_flux.h"

#include <algorithm>
#include <stdexcept>

// ---------------------------------------------------------------------- //
bool FluxRequest::any() const {
     return density || order_parameter || momentum;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::density_x_slot() const {
     return 0;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::density_y_slot() const {
     return 1;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::order_parameter_x_slot(int order_parameter) const {
     check_order_parameter(order_parameter);
     return 2 + 2 * order_parameter;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::order_parameter_y_slot(int order_parameter) const {
     check_order_parameter(order_parameter);
     return 2 + 2 * order_parameter + 1;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::momentum_base_slot() const {
     return 2 + 2 * num_order_parameters_;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::momentum_xx_slot() const {
     return momentum_base_slot();
}
// ---------------------------------------------------------------------- //
int FluxBuffer::momentum_xy_slot() const {
     return momentum_base_slot() + 1;
}
// ---------------------------------------------------------------------- //
// ---------------------------------------------------------------------- //
int FluxBuffer::momentum_yy_slot() const {
     return momentum_base_slot() + 2;
}
// ---------------------------------------------------------------------- //
int FluxBuffer::num_slots() const {
     return 2 + 2 * num_order_parameters_ + 3;
}
// ---------------------------------------------------------------------- //
std::size_t FluxBuffer::offset(int slot) const {
     return static_cast<std::size_t>(slot) * local_spectral_size_;
}
// ---------------------------------------------------------------------- //
void FluxBuffer::check_order_parameter(int order_parameter) const {
     if (order_parameter < 0 || order_parameter >= num_order_parameters_) {
          throw std::out_of_range("FluxBuffer order-parameter index out of range.");
     }
}
// ---------------------------------------------------------------------- //
void FluxBuffer::check_compatible(const FluxBuffer& flux) const {
     if (num_order_parameters_ != flux.num_order_parameters_
         || local_spectral_size_ != flux.local_spectral_size_) {
          throw std::runtime_error("FluxBuffer size mismatch.");
     }
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::slot_data(int slot) const {
     return data_.data() + offset(slot);
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::slot_data(int slot) {
     return data_.data() + offset(slot);
}
// ---------------------------------------------------------------------- //
void FluxBuffer::accumulate_slot(int slot, double weight, const Complex* values) {
     if (weight == 0.0) {
          return;
     }

     Complex* dst = slot_data(slot);
     for (std::size_t i = 0; i < local_spectral_size_; ++i) {
          dst[i] += weight * values[i];
     }
}
// ---------------------------------------------------------------------- //
FluxBuffer::FluxBuffer(const Params& params, const Domain2D& domain)
     : num_order_parameters_(params.physics.num_order_parameters),
       local_spectral_size_(domain.spectral_size()) {
     if (num_order_parameters_ < 0) {
          throw std::runtime_error("FluxBuffer requires a nonnegative number of order parameters.");
     }

     data_.assign(static_cast<std::size_t>(num_slots()) * local_spectral_size_, Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
void FluxBuffer::set_request(const FluxRequest& request) {
     request_ = request;
}
// ---------------------------------------------------------------------- //
const FluxRequest& FluxBuffer::request() const {
     return request_;
}
// ---------------------------------------------------------------------- //
bool FluxBuffer::requested() const {
     return request_.any();
}
// ---------------------------------------------------------------------- //
void FluxBuffer::begin_step() {
     if (request_.any()) {
          std::fill(data_.begin(), data_.end(), Complex(0.0, 0.0));
     }
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::density_flux_x_hat_data() {
     return slot_data(density_x_slot());
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::density_flux_x_hat_data() const {
     return slot_data(density_x_slot());
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::density_flux_y_hat_data() {
     return slot_data(density_y_slot());
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::density_flux_y_hat_data() const {
     return slot_data(density_y_slot());
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::order_parameter_flux_x_hat_data(int order_parameter) {
     return slot_data(order_parameter_x_slot(order_parameter));
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::order_parameter_flux_x_hat_data(int order_parameter) const {
     return slot_data(order_parameter_x_slot(order_parameter));
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::order_parameter_flux_y_hat_data(int order_parameter) {
     return slot_data(order_parameter_y_slot(order_parameter));
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::order_parameter_flux_y_hat_data(int order_parameter) const {
     return slot_data(order_parameter_y_slot(order_parameter));
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::momentum_flux_xx_hat_data() {
     return slot_data(momentum_xx_slot());
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::momentum_flux_xx_hat_data() const {
     return slot_data(momentum_xx_slot());
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::momentum_flux_xy_hat_data() {
     return slot_data(momentum_xy_slot());
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::momentum_flux_xy_hat_data() const {
     return slot_data(momentum_xy_slot());
}
// ---------------------------------------------------------------------- //
Complex* FluxBuffer::momentum_flux_yy_hat_data() {
     return slot_data(momentum_yy_slot());
}
// ---------------------------------------------------------------------- //
const Complex* FluxBuffer::momentum_flux_yy_hat_data() const {
     return slot_data(momentum_yy_slot());
}
// ---------------------------------------------------------------------- //
void FluxBuffer::accumulate_stage_flux(double weight, const FluxBuffer& stage_flux) {
     if (!request_.any() || weight == 0.0) {
          return;
     }

     check_compatible(stage_flux);

     if (request_.density) {
          accumulate_slot(density_x_slot(), weight, stage_flux.density_flux_x_hat_data());
          accumulate_slot(density_y_slot(), weight, stage_flux.density_flux_y_hat_data());
     }

     if (request_.order_parameter) {
          for (int op = 0; op < num_order_parameters_; ++op) {
               accumulate_slot(order_parameter_x_slot(op), weight, stage_flux.order_parameter_flux_x_hat_data(op));
               accumulate_slot(order_parameter_y_slot(op), weight, stage_flux.order_parameter_flux_y_hat_data(op));
          }
     }

     if (request_.momentum) {
          accumulate_slot(momentum_xx_slot(), weight, stage_flux.momentum_flux_xx_hat_data());
          accumulate_slot(momentum_xy_slot(), weight, stage_flux.momentum_flux_xy_hat_data());
          accumulate_slot(momentum_yy_slot(), weight, stage_flux.momentum_flux_yy_hat_data());
     }
}
// ---------------------------------------------------------------------- //
void FluxBuffer::add_incompressible_projection_pressure_flux(const SpectralMask2D& spectral_mask) {
    if (!request_.momentum) {
        return;
    }

    Complex* pi_xx = momentum_flux_xx_hat_data();
    Complex* pi_xy = momentum_flux_xy_hat_data();
    Complex* pi_yy = momentum_flux_yy_hat_data();

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.k2 == 0.0) {
            continue;
        }

        const std::size_t i = mode.index;
        const double kx = mode.kx;
        const double ky = mode.ky;

        const Complex rhs_x =
            -Complex(0.0, 1.0) * (kx * pi_xx[i] + ky * pi_xy[i]);
        const Complex rhs_y =
            -Complex(0.0, 1.0) * (kx * pi_xy[i] + ky * pi_yy[i]);

        const Complex k_dot_rhs = kx * rhs_x + ky * rhs_y;
        const Complex pressure_hat = -Complex(0.0, 1.0) * k_dot_rhs / mode.k2;

        pi_xx[i] += pressure_hat;
        pi_yy[i] += pressure_hat;
    }
}
// ---------------------------------------------------------------------- //
