CXX := mpic++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic -Isrc -I/usr/local/include
LDFLAGS := -L/usr/local/lib
LDLIBS := -lheffte \
          -lfftw3_mpi -lfftw3 \
          -lfftw3f_mpi -lfftw3f \
          -lfftw3l_mpi -lfftw3l

TARGET := out.exe
SRCS := \
	src/main.cc \
	src/buffer_physical_state.cc \
	src/domain.cc \
	src/initial_condition_registry_builtin.cc \
	src/initial_condition_uniform_density.cc \
	src/initial_condition_uniform_momentum.cc \
	src/measure_registry_builtin.cc \
	src/measure_snapshot.cc \
	src/monitor.cc \
	src/param_parser.cc \
	src/state.cc \
	src/simulationinfo.cc

INPUT ?= examples/01_two_component_ideal_gas/input.script

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	./$(TARGET) $(INPUT)

clean:
	rm -f $(TARGET)
