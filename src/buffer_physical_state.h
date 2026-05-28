#ifndef SFI_BUFFER_PHYSICAL_STATE_H
#define SFI_BUFFER_PHYSICAL_STATE_H

#include "domain.h"
#include "simulationinfo.h"

#include <cstddef>
#include <vector>

class PhysicalStateBuffer {
private:
    const Domain2D& domain_;
    int num_order_parameters_ = 0;
    int num_fields_ = 0;
    std::size_t local_physical_size_ = 0;
    std::vector<double> physical_data_;

    std::size_t field_offset(int field_index) const;
    std::size_t rho_offset() const { return 0; }
    std::size_t psi_offset(int order_parameter) const;
    std::size_t jx_offset() const;
    std::size_t jy_offset() const;

    void check_order_parameter(int order_parameter) const;

public:
    PhysicalStateBuffer(const Params& params, const Domain2D& domain);

    double& rho(int gx, int gy);
    const double& rho(int gx, int gy) const;

    double& psi(int order_parameter, int gx, int gy);
    const double& psi(int order_parameter, int gx, int gy) const;

    double& jx(int gx, int gy);
    const double& jx(int gx, int gy) const;

    double& jy(int gx, int gy);
    const double& jy(int gx, int gy) const;

    double* rho_data();
    const double* rho_data() const;
    double* psi_data(int order_parameter);
    const double* psi_data(int order_parameter) const;
    double* jx_data();
    const double* jx_data() const;
    double* jy_data();
    const double* jy_data() const;

    double* data();
    const double* data() const;
};

#endif
