#ifndef SHHD_MEASURE_FLUX_H
#define SHHD_MEASURE_FLUX_H

#include "measure.h"

#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class FluxMeasure : public Measure {
private:
    struct FluxField {
        std::string name;
        const Complex* data = nullptr;
    };

    int num_order_parameters_ = 0;
    int nevery_ = 0;
    std::string file_;
    FluxRequest request_;

    std::ofstream output_;
    void open_output_if_needed();

    std::vector<std::string> selected_field_names() const;
    std::vector<FluxField> selected_fields(const FluxBuffer& flux) const;
    std::optional<std::size_t> zero_mode_index() const;

public:
    FluxMeasure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command);

    FluxRequest flux_request() const override;
    void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
    void finalize() override;
};

#endif
