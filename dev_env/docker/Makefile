CXX := mpic++
CXXFLAGS := -std=c++14 -O2 -Wall -Wextra -pedantic -I/usr/include
LDFLAGS := -L/lib/aarch64-linux-gnu
LDLIBS := -lfftw3_mpi -lfftw3

TARGET := out.exe
SRC := check_mpi_fftw.cpp

.PHONY: all run paths clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	mpirun --allow-run-as-root -np 4 ./$(TARGET)

paths:
	@echo "Compiler wrapper:"
	@which $(CXX)
	@echo
	@echo "MPI compile flags:"
	@$(CXX) --showme:compile
	@echo
	@echo "MPI link flags:"
	@$(CXX) --showme:link
	@echo
	@echo "FFTW headers:"
	@ls -l /usr/include/fftw3.h /usr/include/fftw3-mpi.h
	@echo
	@echo "FFTW libraries:"
	@ldconfig -p | grep 'libfftw3'

clean:
	rm -f $(TARGET)
