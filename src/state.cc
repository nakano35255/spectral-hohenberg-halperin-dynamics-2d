#include "state.h"

#include <algorithm>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
State::State(const Domain2D& domain, const Params& params)
    : domain_(domain),
      num_order_parameters_(params.physics.num_order_parameters),
      num_fields_(params.physics.num_order_parameters + 3),
      local_spectral_size_(domain.spectral_size()) {
    if (num_order_parameters_ < 0) {
        throw std::runtime_error("State requires a nonnegative number of order parameters.");
    }

    spectral_data_.assign(
        static_cast<std::size_t>(num_fields_) * local_spectral_size_,
        Complex(0.0, 0.0)
    );
}
// ---------------------------------------------------------------------- //
std::size_t State::field_offset(int field_index) const {
    return static_cast<std::size_t>(field_index) * local_spectral_size_;
}
// ---------------------------------------------------------------------- //
std::size_t State::psi_offset(int order_parameter) const {
    check_order_parameter(order_parameter);
    return field_offset(order_parameter + 1);
}
// ---------------------------------------------------------------------- //
std::size_t State::jx_offset() const {
    return field_offset(num_order_parameters_ + 1);
}
// ---------------------------------------------------------------------- //
std::size_t State::jy_offset() const {
    return field_offset(num_order_parameters_ + 2);
}
// ---------------------------------------------------------------------- //
void State::check_order_parameter(int order_parameter) const {
    if (order_parameter < 0 || order_parameter >= num_order_parameters_) {
        throw std::out_of_range("State order-parameter index out of range: " + std::to_string(order_parameter));
    }
}
// ---------------------------------------------------------------------- //
Complex& State::rho_hat(int gx, int gy) {
    return spectral_data_[rho_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::rho_hat(int gx, int gy) const {
    return spectral_data_[rho_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex& State::psi_hat(int order_parameter, int gx, int gy) {
    return spectral_data_[psi_offset(order_parameter) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::psi_hat(int order_parameter, int gx, int gy) const {
    return spectral_data_[psi_offset(order_parameter) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex& State::jx_hat(int gx, int gy) {
    return spectral_data_[jx_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::jx_hat(int gx, int gy) const {
    return spectral_data_[jx_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex& State::jy_hat(int gx, int gy) {
    return spectral_data_[jy_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::jy_hat(int gx, int gy) const {
    return spectral_data_[jy_offset() + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex* State::rho_hat_data() {
    return spectral_data_.data() + rho_offset();
}
// ---------------------------------------------------------------------- //
const Complex* State::rho_hat_data() const {
    return spectral_data_.data() + rho_offset();
}
// ---------------------------------------------------------------------- //
Complex* State::psi_hat_data(int order_parameter) {
    return spectral_data_.data() + psi_offset(order_parameter);
}
// ---------------------------------------------------------------------- //
const Complex* State::psi_hat_data(int order_parameter) const {
    return spectral_data_.data() + psi_offset(order_parameter);
}
// ---------------------------------------------------------------------- //
Complex* State::jx_hat_data() {
    return spectral_data_.data() + jx_offset();
}
// ---------------------------------------------------------------------- //
const Complex* State::jx_hat_data() const {
    return spectral_data_.data() + jx_offset();
}
// ---------------------------------------------------------------------- //
Complex* State::jy_hat_data() {
    return spectral_data_.data() + jy_offset();
}
// ---------------------------------------------------------------------- //
const Complex* State::jy_hat_data() const {
    return spectral_data_.data() + jy_offset();
}
// ---------------------------------------------------------------------- //
Complex* State::data() {
    return spectral_data_.data();
}
// ---------------------------------------------------------------------- //
const Complex* State::data() const {
    return spectral_data_.data();
}
// ---------------------------------------------------------------------- //
void State::clear() {
    std::fill(spectral_data_.begin(), spectral_data_.end(), Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
