CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic -Isrc
LDFLAGS :=
LDLIBS :=

TARGET := out.exe
SRCS := \
	src/main.cc \
	src/param_parser.cc \
	src/simulationinfo.cc
OBJS := $(SRCS:.cc=.o)

INPUT ?= examples/01_two_component_ideal_gas/input.script

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) $(INPUT)

clean:
	rm -f $(TARGET) $(OBJS)
