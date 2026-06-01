#ifndef SHHD_MEASURE_YOKOTA_GREEN_KUBO_H
#define SHHD_MEASURE_YOKOTA_GREEN_KUBO_H

#include "measure.h"
#include "spectral_mask.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

class YokotaGreenKuboMeasure : public Measure {
private:
    int nevery_ = 0;
    int nblock_ = 0;
    int nsample_ = 0;
    int ntau_ = 0;
    int max_nx_ = 0;
    int max_ny_ = 0;

    std::string file_;
    std::string mode_selection_;
    double kBT_ = 1.0;
    int block_step_ = 0;
    std::size_t nmodes_ = 0;

    SpectralMask2D spectral_mask_;

    std::pair<int, int> wave_number(std::size_t mode_index) const;
    std::vector<std::optional<std::size_t>> owned_local_indices_;
    std::vector<Complex> integrals_;
    std::vector<double> sum_by_mode_tau_;

    int wrap_gy(int gy) const;
    std::optional<std::size_t> local_mode_index(int nx, int ny) const;
    std::size_t mode_tau_index(std::size_t mode_index, int tau_index) const;
    bool write_output_table() const;    

public:
    YokotaGreenKuboMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command
    );

    FluxRequest flux_request() const override;
    void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
};

#endif
