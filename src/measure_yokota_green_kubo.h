#ifndef SHHD_MEASURE_YOKOTA_GREEN_KUBO_H
#define SHHD_MEASURE_YOKOTA_GREEN_KUBO_H

#include "measure.h"
#include "spectral_mask.h"

#include <complex>
#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

class YokotaGreenKuboMeasure : public Measure {
private:
    struct ModeAccumulator {
        int requested_nx = 0;
        int requested_ny = 0;
        int stored_gx = 0;
        int stored_gy = 0;
        bool use_conjugate = false;
        bool owns_mode = false;
        std::size_t mode_index = 0;
        double kx = 0.0;
        double ky = 0.0;
        Complex integral = Complex(0.0, 0.0);
        std::vector<double> sum_by_tau;
        std::vector<int> count_by_tau;
    };

    int nevery_ = 0;
    int nblock_ = 0;
    std::string file_;
    std::string mode_selection_;
    double kBT_ = 1.0;

    SpectralMask2D spectral_mask_;

    int block_step_ = 0;
    int block_index_ = 0;
    std::vector<ModeAccumulator> modes_;

    std::ofstream output_;

    void configure_modes();
    void add_mode(int nx, int ny);
    void open_output_if_needed();
    int wrap_gy(int gy) const;
    double cell_area() const;
    double system_area() const;

public:
    YokotaGreenKuboMeasure(
        const Params& params,
        const Domain2D& domain,
        std::shared_ptr<const MeasureCommandBase> command
    );

    FluxRequest flux_request() const override;
    void observe(
        const State& state,
        FourierTransform2D& fft,
        MeasureWorkspace& workspace,
        const FluxBuffer& flux,
        int step,
        double time
    ) override;
    void finalize() override;
};

#endif
