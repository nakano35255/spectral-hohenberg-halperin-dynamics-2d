#include "buffer_physical_state.h"

#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
PhysicalStateBuffer::PhysicalStateBuffer(const Domain2D& domain, const Params& params)
    : domain_(domain),
      num_components_(params.physics.num_components),
      local_physical_size_(domain.physical_size()) {
    if (num_components_ <= 0) {
        throw std::runtime_error("PhysicalStateBuffer requires a positive number of components.");
    }

    physical_data_.assign(
        static_cast<std::size_t>(num_components_ + 2) * local_physical_size_, 0.0
    );
}
// ---------------------------------------------------------------------- //
std::size_t PhysicalStateBuffer::field_offset(int field_index) const {
    return static_cast<std::size_t>(field_index) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
void PhysicalStateBuffer::check_component(int component) const {
    if (component < 0 || component >= num_components_) {
        throw std::out_of_range("PhysicalStateBuffer component index out of range: " + std::to_string(component));
    }
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::rho(int component, int gx, int gy) {
    check_component(component);
    return physical_data_[field_offset(component) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::rho(int component, int gx, int gy) const {
    check_component(component);
    return physical_data_[field_offset(component) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::jx(int gx, int gy) {
    return physical_data_[field_offset(num_components_) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::jx(int gx, int gy) const {
    return physical_data_[field_offset(num_components_) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::jy(int gx, int gy) {
    return physical_data_[field_offset(num_components_ + 1) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::jy(int gx, int gy) const {
    return physical_data_[field_offset(num_components_ + 1) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::rho_data(int component) {
    check_component(component);
    return physical_data_.data() + field_offset(component);
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::jx_data() {
    return physical_data_.data() + field_offset(num_components_);
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::jy_data() {
    return physical_data_.data() + field_offset(num_components_ + 1);
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::data() {
    return physical_data_.data();
}
// ---------------------------------------------------------------------- //
const double* PhysicalStateBuffer::data() const {
    return physical_data_.data();
}
// ---------------------------------------------------------------------- //
