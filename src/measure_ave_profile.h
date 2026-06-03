#ifndef SHHD_MEASURE_AVE_PROFILE_H
#define SHHD_MEASURE_AVE_PROFILE_H

#include "measure.h"
#include "fcalculator_dynamics_mode.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

enum class AveProfileTargetKind {
     Rho, Psi, Jx, Jy, Vx, Vy,
     OrderParameterFluxX, OrderParameterFluxY,
     MomentumFluxXX, MomentumFluxXY, MomentumFluxYY
};

enum class AveProfileAverageMode {
     Block,
     Running
};

struct AveProfileTarget {
     AveProfileTargetKind kind;
     int component = -1;
     std::string name;
};

class AveProfileMeasure : public Measure {
private:
     int axis_ = 1;
     int nevery_ = 0;
     int nblock_ = 0;
     int profile_size_ = 0;
     int block_step_ = 0;
     int samples_in_block_ = 0;
     int running_samples_ = 0;
     int completed_blocks_ = 0;

     std::string file_;
     AveProfileAverageMode average_mode_ = AveProfileAverageMode::Block;
     std::vector<AveProfileTarget> targets_;
     DynamicsMode dynamics_mode_;

     // for spectral data
     std::size_t local_spectral_size_ = 0;
     int num_spectral_output_fields_ = 0;
     std::vector<int> spectral_slot_by_target_;
     std::vector<Complex> spectral_block_sum_;
     std::vector<Complex> spectral_running_sum_;

     // for physical data (only velocity)
     std::size_t local_physical_size_ = 0;
     int compressible_vx_slot_ = -1;
     int compressible_vy_slot_ = -1;
     int num_compressible_velocity_fields_ = 0;
     std::vector<double> compressible_velocity_block_sum_;
     std::vector<double> compressible_velocity_running_sum_;

     std::ofstream block_output_;

     bool needs_incompressible_velocity() const;
     bool needs_order_parameter_flux() const;
     bool needs_momentum_flux() const;
     bool is_compressible_velocity_target(const AveProfileTarget& target) const;

     bool reference_density_ready_ = false;
     double reference_density_ = 0.0;
     double reference_density(const State& state);

     void add_to_spectral_block(int slot, const Complex* values, double scale = 1.0);
     void add_to_compressible_velocity_block(int slot, const PhysicalStateBuffer& physical, int component);

     void accumulate_sample(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux);
     void finish_block(FourierTransform2D& fft, int step, double time);

     void open_block_output_if_needed();
     void write_header(std::ostream& out) const;
     void write_rows(std::ostream& out, int step, double time, int samples, const std::vector<double>& values) const;

public:
     AveProfileMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command);

     FluxRequest flux_request() const override;
     void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
     void finalize() override;
};

#endif
