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
    PhysicalStateBuffer(const Domain2D& domain, const Params& params);

    double& rho(int gx, int gy);
    const double& rho(int gx, int gy) const;

    double& psi(int order_parameter, int gx, int gy);
    const double& psi(int order_parameter, int gx, int gy) const;

    double& jx(int gx, int gy);
    const double& jx(int gx, int gy) const;

    double& jy(int gx, int gy);
    const double& jy(int gx, int gy) const;

    double* rho_data();
    double* psi_data(int order_parameter);
    double* jx_data();
    double* jy_data();

    double* data();
    const double* data() const;
};

#endif
