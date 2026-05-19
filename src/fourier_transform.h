#ifndef SFI_FOURIER_TRANSFORM_H
#define SFI_FOURIER_TRANSFORM_H

#include "domain.h"

#include <complex>
#include <cstddef>
#include <stdexcept>

class FourierTransform2D {
public:
     using fft_type = heffte::fft3d_r2c<heffte::backend::fftw>;

private:
     const Domain2D& domain_;
     fft_type fft_;

public:
     explicit FourierTransform2D(const Domain2D& domain)
          : domain_(domain),
               fft_(domain.physical_box(),
                    domain.spectral_box(),
                    Domain2D::R2C_DIRECTION,
                    domain.comm()) {}

     void forward(const double* physical, Complex* spectral) {
          forward_many(1, physical, spectral);
     }

     void backward(const Complex* spectral, double* physical) {
          backward_many(1, spectral, physical);
     }

     void forward_many(int nfields, const double* physical, Complex* spectral) {
          if (nfields <= 0) {
               throw std::runtime_error("FourierTransform2D::forward_many requires positive nfields.");
          }
          fft_.forward(nfields, physical, spectral, heffte::scale::none);
     }

     void backward_many(int nfields, const Complex* spectral, double* physical) {
          if (nfields <= 0) {
               throw std::runtime_error("FourierTransform2D::backward_many requires positive nfields.");
          }
          fft_.backward(nfields, spectral, physical, heffte::scale::full);
     }
};

#endif
