#ifndef SFI_SPECTRAL_MASK_H
#define SFI_SPECTRAL_MASK_H

#include "domain.h"
#include "simulationinfo.h"

#include <complex>
#include <cstddef>
#include <vector>

struct SpectralMode2D {
     std::size_t index = 0;
     int gx = 0;
     int gy = 0;
     double kx = 0.0;
     double ky = 0.0;
     double k2 = 0.0;
};

class SpectralMask2D {
private:
     const Domain2D& domain_;
     int active_nx_ = 0;
     int active_ny_ = 0;
     std::vector<SpectralMode2D> active_modes_;
     std::vector<std::size_t> inactive_indices_;

     int signed_ky(int gy) const;
     bool is_active_mode(int gx, int gy) const;

public:
     SpectralMask2D(const Params& params, const Domain2D& domain);

     const std::vector<SpectralMode2D>& active_modes() const {
          return active_modes_;
     }

     bool active(int gx, int gy) const {
          return is_active_mode(gx, gy);
     }

     void apply(Complex* field) const;
     void apply_many(int nfields, Complex* fields) const;
};

#endif
