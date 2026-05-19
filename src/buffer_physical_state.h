#ifndef SFI_BUFFER_PHYSICAL_STATE_H
#define SFI_BUFFER_PHYSICAL_STATE_H

#include "domain.h"
#include "simulationinfo.h"

#include <cstddef>
#include <vector>

class PhysicalStateBuffer {
private:
    const Domain2D& domain_;
    int num_components_ = 0;
    std::size_t local_physical_size_ = 0;
    std::vector<double> physical_data_;

    std::size_t field_offset(int field_index) const;
    void check_component(int component) const;

public:
    PhysicalStateBuffer(const Domain2D& domain, const Params& params);

    double& rho(int component, int gx, int gy);
    const double& rho(int component, int gx, int gy) const;

    double& jx(int gx, int gy);
    const double& jx(int gx, int gy) const;

    double& jy(int gx, int gy);
    const double& jy(int gx, int gy) const;

    double* rho_data(int component);
    double* jx_data();
    double* jy_data();

    double* data();
    const double* data() const;
};

#endif
