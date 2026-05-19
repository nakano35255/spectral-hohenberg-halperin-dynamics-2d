#include "state.h"

#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
State::State(const Domain2D& domain, const Params& params)
    : domain_(domain),
      num_components_(params.physics.num_components),
      local_spectral_size_(domain.spectral_size()) {
    if (num_components_ <= 0) {
        throw std::runtime_error("State requires a positive number of components.");
    }

    spectral_data_.assign(
        static_cast<std::size_t>(num_components_ + 2) * local_spectral_size_, Complex(0.0, 0.0)
    );
}
// ---------------------------------------------------------------------- //
std::size_t State::field_offset(int field_index) const {
    return static_cast<std::size_t>(field_index) * local_spectral_size_;
}
// ---------------------------------------------------------------------- //
void State::check_component(int component) const {
    if (component < 0 || component >= num_components_) {
        throw std::out_of_range("State component index out of range: " + std::to_string(component));
    }
}
// ---------------------------------------------------------------------- //
Complex& State::rho_hat(int component, int gx, int gy) {
    check_component(component);
    return spectral_data_[field_offset(component) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::rho_hat(int component, int gx, int gy) const {
    check_component(component);
    return spectral_data_[field_offset(component) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex& State::jx_hat(int gx, int gy) {
    return spectral_data_[field_offset(num_components_) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::jx_hat(int gx, int gy) const {
    return spectral_data_[field_offset(num_components_) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex& State::jy_hat(int gx, int gy) {
    return spectral_data_[field_offset(num_components_ + 1) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const Complex& State::jy_hat(int gx, int gy) const {
    return spectral_data_[field_offset(num_components_ + 1) + domain_.spectral_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
Complex* State::rho_hat_data(int component) {
    check_component(component);
    return spectral_data_.data() + field_offset(component);
}
// ---------------------------------------------------------------------- //
Complex* State::jx_hat_data() {
    return spectral_data_.data() + field_offset(num_components_);
}
// ---------------------------------------------------------------------- //
Complex* State::jy_hat_data() {
    return spectral_data_.data() + field_offset(num_components_ + 1);
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
