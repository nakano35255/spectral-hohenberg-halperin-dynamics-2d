#include "spectral_mask.h"

#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------------- //
SpectralMask2D::SpectralMask2D(const Domain2D& domain, const Params& params)
     : domain_(domain),
       active_nx_(params.grid.active_num_grid[0]),
       active_ny_(params.grid.active_num_grid[1])
{
     if (active_nx_ <= 0 || active_ny_ <= 0) {
          throw std::runtime_error("SpectralMask2D requires positive active grid sizes.");
     }
     if (active_nx_ > domain_.nx_global() || active_ny_ > domain_.ny_global()) {
          throw std::runtime_error("SpectralMask2D active grid cannot exceed compute grid.");
     }

     const Box2D& box = domain_.spectral_box();

     std::size_t index = 0;
     for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
          const double ky = domain_.ky(gy);

          for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
               const double kx = domain_.kx(gx);
               const double k2 = kx * kx + ky * ky;

               if (is_active_mode(gx, gy)) {
                    active_modes_.push_back({index, gx, gy, kx, ky, k2});
               } else {
                    inactive_indices_.push_back(index);
               }

               ++index;
          }
     }
}
// ---------------------------------------------------------------------- //
int SpectralMask2D::signed_ky(int gy) const {
     const int ny = domain_.ny_global();
     return (gy <= ny / 2) ? gy : gy - ny;
}
// ---------------------------------------------------------------------- //
bool SpectralMask2D::is_active_mode(int gx, int gy) const {
     if (gx < 0 || gx > active_nx_ / 2) {
          return false;
     }

     const int ky = signed_ky(gy);
     return std::abs(ky) <= active_ny_ / 2;
}
// ---------------------------------------------------------------------- //
void SpectralMask2D::apply(Complex* field) const {
     for (std::size_t index : inactive_indices_) {
          field[index] = Complex(0.0, 0.0);
     }
}
// ---------------------------------------------------------------------- //
void SpectralMask2D::apply_many(int nfields, Complex* fields) const {
     if (nfields <= 0) {
          throw std::runtime_error("SpectralMask2D::apply_many requires positive nfields.");
     }

     const std::size_t local_spectral_size = domain_.spectral_size();

     for (int field = 0; field < nfields; ++field) {
          Complex* data = fields + static_cast<std::size_t>(field) * local_spectral_size;
          apply(data);
     }
}
// ---------------------------------------------------------------------- //
