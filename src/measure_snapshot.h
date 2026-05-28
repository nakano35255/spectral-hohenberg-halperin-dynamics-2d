#ifndef MEASURE_SNAPSHOT_H
#define MEASURE_SNAPSHOT_H

#include <memory>
#include <string>
#include <vector>

#include "measure.h"

class SnapshotMeasure: public Measure {
private:
     int num_order_parameters_ = 0;
     int num_fields_ = 0;
     int nevery_ = 0;
     std::string file_;
     std::string space_;
     int step_width_ = 1;

     static int compute_step_width(int max_steps);
     static std::string measurement_step_string(int value, int width);
     std::string snapshot_filename(int step) const;
     void write_physical_snapshot(const std::string& filename, const std::vector<double>& global_data, const Domain2D& domain, int step, double time) const;
     void write_spectral_snapshot(const std::string& filename, const std::vector<double>& global_data, const Domain2D& domain, int step, double time) const;
     void observe_physical(PhysicalStateBuffer& physical, const Domain2D& domain, int step, double time) const;
     void observe_spectral(const State& state, const Domain2D& domain, int step, double time) const;
    
public:
     SnapshotMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command);
     void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, int step, double time) override;
};

#endif
