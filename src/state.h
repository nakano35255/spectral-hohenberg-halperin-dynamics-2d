#ifndef SFI_STATE_H
#define SFI_STATE_H

#include "domain.h"
#include "simulationinfo.h"

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

class State {
private:
    const Domain2D& domain_;
    int num_order_parameters_ = 0;
    int num_fields_ = 0;
    std::size_t local_spectral_size_ = 0;
    std::vector<Complex> spectral_data_;

    std::size_t field_offset(int field_index) const;

    std::size_t rho_offset() const { return 0; }
    std::size_t psi_offset(int order_parameter) const;
    std::size_t jx_offset() const;
    std::size_t jy_offset() const;

    void check_order_parameter(int order_parameter) const;

public:
    State(const Params& params, const Domain2D& domain);

    int num_order_parameters() const { return num_order_parameters_; }
    int num_fields() const { return num_fields_; }

    Complex& rho_hat(int gx, int gy);
    const Complex& rho_hat(int gx, int gy) const;

    Complex& psi_hat(int order_parameter, int gx, int gy);
    const Complex& psi_hat(int order_parameter, int gx, int gy) const;

    Complex& jx_hat(int gx, int gy);
    const Complex& jx_hat(int gx, int gy) const;

    Complex& jy_hat(int gx, int gy);
    const Complex& jy_hat(int gx, int gy) const;

    Complex* rho_hat_data();
    const Complex* rho_hat_data() const;

    Complex* psi_hat_data(int order_parameter);
    const Complex* psi_hat_data(int order_parameter) const;

    Complex* jx_hat_data();
    const Complex* jx_hat_data() const;

    Complex* jy_hat_data();
    const Complex* jy_hat_data() const;

    Complex* data();
    const Complex* data() const;

    void clear();
};

#endif
