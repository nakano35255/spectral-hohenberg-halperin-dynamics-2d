#ifndef SHHD_MEASURE_ENERGETICS_H
#define SHHD_MEASURE_ENERGETICS_H

#include "fcalculator_dynamics_mode.h"
#include "measure.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "spectral_mask.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

class EnergeticsMeasure : public Measure {
private:
    const Thermodynamics& thermodynamics_;
    const FreeEnergy& free_energy_;
    const TransportCoefficient& transport_;
    DynamicsMode dynamics_mode_;
    SpectralMask2D spectral_mask_;

    int num_order_parameters_ = 0;
    int num_fields_ = 0;
    int nevery_ = 0;
    std::string file_;

    std::vector<double> k0_coefficients_;
    std::vector<double> k2_coefficients_;
    std::vector<double> k4_coefficients_;
    double linear_pressure_coefficient_ = 0.0;

    std::ofstream output_;
    void open_output_if_needed();

    double compute_local_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_free_energy(const State& state, const PhysicalStateBuffer* physical) const;
    // helper
    double cell_area() const;
    double reference_density(const State& state) const;
    double compute_local_compressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_incompressible_kinetic_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_quadratic_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_physical_compressibility_energy(const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_quadratic_free_energy(int order_parameter, const State& state, const PhysicalStateBuffer* physical) const;
    double compute_local_physical_free_energy(const State& state, const PhysicalStateBuffer* physical) const;

public:
    EnergeticsMeasure(const Params& params, const Domain2D& domain, const Thermodynamics& thermodynamics, const FreeEnergy& free_energy, const TransportCoefficient& transport, std::shared_ptr<const MeasureCommandBase> command);
    void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) override;

    void finalize() override;
};

#endif
