#include "buffer_physical_state.h"

#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
PhysicalStateBuffer::PhysicalStateBuffer(const Params& params, const Domain2D& domain)
    : domain_(domain),
      num_order_parameters_(params.physics.num_order_parameters),
      num_fields_(params.physics.num_order_parameters + 3),
      local_physical_size_(domain.physical_size()) {
    if (num_order_parameters_ < 0) {
        throw std::runtime_error("PhysicalStateBuffer requires a nonnegative number of order parameters.");
    }

    physical_data_.assign(
        static_cast<std::size_t>(num_fields_) * local_physical_size_, 0.0
    );
}
// ---------------------------------------------------------------------- //
std::size_t PhysicalStateBuffer::field_offset(int field_index) const {
    return static_cast<std::size_t>(field_index) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
std::size_t PhysicalStateBuffer::psi_offset(int order_parameter) const {
    check_order_parameter(order_parameter);
    return field_offset(order_parameter + 1);
}
// ---------------------------------------------------------------------- //
std::size_t PhysicalStateBuffer::jx_offset() const {
    return field_offset(num_order_parameters_ + 1);
}
// ---------------------------------------------------------------------- //
std::size_t PhysicalStateBuffer::jy_offset() const {
    return field_offset(num_order_parameters_ + 2);
}
// ---------------------------------------------------------------------- //
void PhysicalStateBuffer::check_order_parameter(int order_parameter) const {
    if (order_parameter < 0 || order_parameter >= num_order_parameters_) {
        throw std::out_of_range("PhysicalStateBuffer order-parameter index out of range: " + std::to_string(order_parameter));
    }
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::rho(int gx, int gy) {
    return physical_data_[rho_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::rho(int gx, int gy) const {
    return physical_data_[rho_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::psi(int order_parameter, int gx, int gy) {
    return physical_data_[psi_offset(order_parameter) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::psi(int order_parameter, int gx, int gy) const {
    return physical_data_[psi_offset(order_parameter) + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::jx(int gx, int gy) {
    return physical_data_[jx_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::jx(int gx, int gy) const {
    return physical_data_[jx_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double& PhysicalStateBuffer::jy(int gx, int gy) {
    return physical_data_[jy_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
const double& PhysicalStateBuffer::jy(int gx, int gy) const {
    return physical_data_[jy_offset() + domain_.physical_local_index(gx, gy)];
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::rho_data() {
    return physical_data_.data() + rho_offset();
}
// ---------------------------------------------------------------------- //
const double* PhysicalStateBuffer::rho_data() const {
    return physical_data_.data() + rho_offset();
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::psi_data(int order_parameter) {
    return physical_data_.data() + psi_offset(order_parameter);
}
// ---------------------------------------------------------------------- //
const double* PhysicalStateBuffer::psi_data(int order_parameter) const {
    return physical_data_.data() + psi_offset(order_parameter);
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::jx_data() {
    return physical_data_.data() + jx_offset();
}
// ---------------------------------------------------------------------- //
const double* PhysicalStateBuffer::jx_data() const {
    return physical_data_.data() + jx_offset();
}
// ---------------------------------------------------------------------- //
double* PhysicalStateBuffer::jy_data() {
    return physical_data_.data() + jy_offset();
}
// ---------------------------------------------------------------------- //
const double* PhysicalStateBuffer::jy_data() const {
    return physical_data_.data() + jy_offset();
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
