CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic -Isrc
LDFLAGS :=
LDLIBS :=

TARGET := out.exe
SRCS := \
	src/main.cc \
	src/measure_registry_builtin.cc \
	src/monitor.cc \
	src/param_parser.cc \
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
