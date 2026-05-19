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
    int num_components_ = 0;
    std::size_t local_spectral_size_ = 0;
    std::vector<Complex> spectral_data_;

    std::size_t field_offset(int field_index) const;

    void check_component(int component) const;

public:
    State(const Domain2D& domain, const Params& params);

    Complex& rho_hat(int component, int gx, int gy);
    const Complex& rho_hat(int component, int gx, int gy) const;

    Complex& jx_hat(int gx, int gy);
    const Complex& jx_hat(int gx, int gy) const;

    Complex& jy_hat(int gx, int gy);
    const Complex& jy_hat(int gx, int gy) const;

    Complex* rho_hat_data(int component);
    Complex* jx_hat_data();
    Complex* jy_hat_data();

    Complex* data();
    const Complex* data() const;

    void clear();
};

#endif
