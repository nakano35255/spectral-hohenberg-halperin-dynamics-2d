#ifndef SHHD_MEASURE_CORRELATION_STATIC_H
#define SHHD_MEASURE_CORRELATION_STATIC_H

#include "measure.h"
#include "spectral_mask.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

enum class CorrelationStaticTargetKind {Rho, Jx, Jy, Psi};
enum class CorrelationStaticOutputMode {TwoD, Shell};
enum class CorrelationStaticAverageMode {Block, Running};

struct CorrelationStaticTarget {
     CorrelationStaticTargetKind kind;
     int component = -1;
     std::string name;
};

class CorrelationStaticMeasure : public Measure {
private:
     int nevery_ = 0;
     int nblock_ = 0;
     int nfields_ = 0;
     int npairs_ = 0;

     int block_step_ = 0;
     int samples_in_block_ = 0;
     int completed_blocks_ = 0;
     int running_samples_ = 0;

     std::string file_;
     CorrelationStaticAverageMode average_mode_ = CorrelationStaticAverageMode::Running;
     CorrelationStaticOutputMode output_mode_ = CorrelationStaticOutputMode::Shell;
     std::vector<CorrelationStaticTarget> targets_;
     bool cross_ = false;

     // for spectral modes
     SpectralMask2D spectral_mask_;
     std::size_t local_spectral_size_ = 0;

     // for shell mode
     double dk_ = 0.0;
     int nshells_ = 0;
     std::vector<double> shell_weight_counts_;
     int shell_index(double k2) const;
     void initialize_shells();

     std::vector<Complex> mode_block_sum_;
     std::vector<Complex> mode_running_sum_;
     std::vector<Complex> shell_block_sum_;
     std::vector<Complex> shell_running_sum_;

     std::size_t pair_mode_index(int pair_index, std::size_t mode_index) const;
     std::size_t pair_shell_index(int pair_index, int shell) const;

     const Complex* target_data(const State& state, const CorrelationStaticTarget& target) const;
     void accumulate_sample(const State& state);

     void finish_block(int step, double time);

     // for output
     std::ofstream block_output_;

     void open_block_output_if_needed();
     void write_header(std::ostream& out) const;
     void write_pair_column_names(std::ostream& out, bool complex_columns) const;
     void write_2d_output(int step, double time, int output_samples, const std::vector<Complex>& source);
     void write_shell_output(int step, double time, int output_samples, const std::vector<Complex>& source);

public:
     CorrelationStaticMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command);
     void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
     void finalize() override;
};

#endif