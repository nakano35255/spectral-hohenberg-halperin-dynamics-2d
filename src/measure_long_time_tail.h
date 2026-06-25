#ifndef SHHD_MEASURE_LONG_TIME_TAIL_H
#define SHHD_MEASURE_LONG_TIME_TAIL_H

#include "measure.h"
#include "spectral_mask.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

enum class LongTimeTailTargetKind {
     Rho, Jx, Jy, Psi
};

enum class LongTimeTailAverageMode {Block, Running};

struct LongTimeTailTarget {
     LongTimeTailTargetKind kind;
     int component = -1;
     std::string name;
};

class LongTimeTailMeasure : public Measure {
private:
     int nevery_ = 0;
     int nblock_ = 0;
     int nfields_ = 0;
     int nlags_ = 0;
     int npairs_ = 0;

     int block_step_ = 0;
     int completed_blocks_ = 0;

     std::string file_;
     LongTimeTailAverageMode average_mode_ = LongTimeTailAverageMode::Running;
     std::vector<LongTimeTailTarget> targets_;
     bool cross_ = false;

     SpectralMask2D spectral_mask_;
     std::size_t local_spectral_size_ = 0;

     std::vector<Complex> block_fields_;
     std::vector<double> block_sum_;
     std::vector<double> running_sum_;

     std::ofstream block_output_;

     std::size_t sample_field_offset(int sample_index, int field_index) const;
     std::size_t pair_lag_index(int pair_index, int lag_index) const;

     void store_current_sample(const State& state, int sample_index);

     double compute_local_correlation(const Complex* later, const Complex* earlier) const;
     void accumulate_auto_correlations(std::vector<double>& local_sum) const;
     void accumulate_cross_correlations(std::vector<double>& local_sum) const;

     void finish_block(int step, double time);
     void reset_block();

     void open_block_output_if_needed();
     void write_header(std::ostream& out) const;
     void write_rows(std::ostream& out, int step, double time, const std::vector<double>& sums) const;

public:
     LongTimeTailMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command);

     void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
     void finalize() override;
};

#endif