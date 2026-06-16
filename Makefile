CXX := mpic++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -pedantic -Isrc -I/usr/local/include
LDFLAGS := -L/usr/local/lib
LDLIBS := -lheffte \
          -lfftw3_mpi -lfftw3 \
          -lfftw3f_mpi -lfftw3f \
          -lfftw3l_mpi -lfftw3l

TARGET := src/out.exe
.DEFAULT_GOAL := all

SRCS := \
	src/main.cc \
	src/buffer_flux.cc \
	src/buffer_physical_state.cc \
	src/domain.cc \
	src/fcalculator.cc \
	src/fcalculator_workspace.cc \
	src/initial_condition.cc \
	src/initial_condition_equilibrium_gaussian_density.cc \
	src/initial_condition_equilibrium_gaussian_momentum_compressible.cc \
	src/initial_condition_equilibrium_gaussian_momentum_incompressible.cc \
	src/initial_condition_equilibrium_gaussian_order_parameter.cc \
	src/initial_condition_registry_builtin.cc \
	src/initial_condition_gaussian_density.cc \
	src/initial_condition_gaussian_order_parameter.cc \
	src/initial_condition_random_vorticity_momentum.cc \
	src/initial_condition_sine_density.cc \
	src/initial_condition_sine_momentum.cc \
	src/initial_condition_sine_order_parameter.cc \
	src/initial_condition_uniform_density.cc \
	src/initial_condition_uniform_momentum.cc \
	src/initial_condition_uniform_order_parameter.cc \
	src/measure_registry_builtin.cc \
	src/measure_manager.cc \
	src/measure_yokota_green_kubo.cc \
	src/measure_snapshot.cc \
	src/measure_time_series.cc \
	src/measure_ave_profile.cc \
	src/model_free_energy_null.cc \
	src/model_free_energy_phi4.cc \
	src/model_free_energy_quadratic.cc \
	src/model_free_energy_registry_builtin.cc \
	src/model_thermodynamics_linear_eos.cc \
	src/model_thermodynamics_null.cc \
	src/model_thermodynamics_registry_builtin.cc \
	src/model_transport_coefficient_constant.cc \
	src/model_transport_coefficient_registry_builtin.cc \
	src/monitor.cc \
	src/param_parser.cc \
	src/spectral_mask.cc \
	src/state.cc \
	src/simulationinfo.cc

PKG_CONFIG := src/packages_enabled.mk
-include $(PKG_CONFIG)

ifeq ($(PKG_PASSIVE_SCALAR),1)
include src/PASSIVE_SCALAR/package.mk
endif


.PHONY: yes-PASSIVE-SCALAR no-PASSIVE-SCALAR package-status

yes-PASSIVE-SCALAR:
	@touch $(PKG_CONFIG)
	@grep -v '^PKG_PASSIVE_SCALAR[[:space:]]*:=' $(PKG_CONFIG) > $(PKG_CONFIG).tmp || true
	@printf '%s\n' 'PKG_PASSIVE_SCALAR := 1' >> $(PKG_CONFIG).tmp
	@mv $(PKG_CONFIG).tmp $(PKG_CONFIG)
	@echo "Enabled package: PASSIVE-SCALAR"
	@echo "Run: make clean && make"

no-PASSIVE-SCALAR:
	@touch $(PKG_CONFIG)
	@grep -v '^PKG_PASSIVE_SCALAR[[:space:]]*:=' $(PKG_CONFIG) > $(PKG_CONFIG).tmp || true
	@mv $(PKG_CONFIG).tmp $(PKG_CONFIG)
	@echo "Disabled package: PASSIVE-SCALAR"
	@echo "Run: make clean && make"

package-status:
	@echo "PASSIVE-SCALAR: $(if $(filter 1,$(PKG_PASSIVE_SCALAR)),yes,no)"


.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f $(TARGET)
