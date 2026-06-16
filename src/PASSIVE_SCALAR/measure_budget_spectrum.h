#ifndef SHHD_PACKAGE_PASSIVE_SCALAR_MEASURE_BUDGET_SPECTRUM_H
#define SHHD_PACKAGE_PASSIVE_SCALAR_MEASURE_BUDGET_SPECTRUM_H

#include "fcalculator_dynamics_mode.h"
#include "measure.h"
#include "model_free_energy.h"
#include "model_transport_coefficient.h"
#include "spectral_mask.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

enum class BudgetSpectrumOutputMode { TwoD, Shell };
enum class BudgetSpectrumAverageMode { Block, Running };

class BudgetSpectrumMeasure : public Measure {
private:
    static constexpr int TERM_TRANSFER = 0;
    static constexpr int TERM_DISSIPATION = 1;
    static constexpr int TERM_PRODUCTION = 2;
    static constexpr int NUM_TERMS = 3;

    const FreeEnergy& free_energy_;
    const TransportCoefficient& transport_;
    DynamicsMode dynamics_mode_;
    SpectralMask2D spectral_mask_;

    int component_ = 0;
    int nevery_ = 1000;
    int nblock_ = 1000;
    int block_step_ = 0;
    int samples_in_block_ = 0;
    int completed_blocks_ = 0;
    int running_samples_ = 0;

    std::string file_;
    BudgetSpectrumOutputMode output_mode_ = BudgetSpectrumOutputMode::TwoD;
    BudgetSpectrumAverageMode average_mode_ = BudgetSpectrumAverageMode::Block;

    int num_order_parameters_ = 0;
    std::size_t local_physical_size_ = 0;
    std::size_t local_spectral_size_ = 0;

    double mobility_ = 0.0;
    double k0_coefficient_ = 0.0;
    double k2_coefficient_ = 0.0;
    double k4_coefficient_ = 0.0;
    bool has_physical_chemical_potential_ = false;

    bool has_transfer_ = false;
    bool has_dissipation_ = false;
    bool has_production_ = false;
    bool need_physical_velocity_ = false;
    bool need_velocity_hat_ = false;

    std::vector<SineForceFixConfig> sine_forces_;
    std::vector<GradientForceFixConfig> gradient_forces_;

    // for shell mode
    double dk_ = 0.0;
    int nshells_ = 0;
    std::vector<double> shell_weight_counts_;
    int shell_index(double k2) const;
    void initialize_shells();

    // for velocity
    double reference_density_ = 0.0;
    bool reference_density_ready_ = false;
    double reference_density(const State& state);
    std::vector<double> velocity_physical_;
    std::vector<Complex> velocity_hat_;
    void prepare_velocity(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace);

    // for transfer
    std::vector<double> transfer_flux_physical_;
    std::vector<Complex> transfer_flux_hat_;
    void prepare_transfer(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace);
    Complex transfer_term(const SpectralMode2D& mode) const;

    // for dissipation
    std::vector<double> chemical_potential_physical_;
    std::vector<Complex> chemical_potential_hat_;
    std::vector<double> psi_point_;
    void prepare_dissipation(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace);
    Complex dissipation_term(const SpectralMode2D& mode, const Complex& psi) const;

    // for production
    Complex production_term(const SpectralMode2D& mode) const;

    std::vector<double> mode_block_sum_;
    std::vector<double> mode_running_sum_;
    std::vector<double> shell_block_sum_;
    std::vector<double> shell_running_sum_;

    void accumulate_sample(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace);

    // for output
    std::ofstream block_output_;

    void finish_block(int step, double time);
    void write_2d_output(int step, double time, int output_samples, const std::vector<double>& source);
    void write_shell_output(int step, double time, int output_samples, const std::vector<double>& source);
    void open_block_output_if_needed();
    void write_header(std::ostream& os) const;

public:
    BudgetSpectrumMeasure(
        const Params& params,
        const Domain2D& domain,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport,
        std::shared_ptr<const MeasureCommandBase> command
    );

    void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;
    void finalize() override;
};

#endif
