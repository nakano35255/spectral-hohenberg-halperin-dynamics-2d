CXX := mpic++
CXXFLAGS := -std=c++14 -O2 -Wall -Wextra -pedantic -I/usr/local/include
LDFLAGS := -L/usr/local/lib
LDLIBS := -lheffte \
          -lfftw3_mpi -lfftw3 \
          -lfftw3f_mpi -lfftw3f \
          -lfftw3l_mpi -lfftw3l

TARGET := out.exe
SRC := check_heffte_fftw.cpp

.PHONY: all run paths clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	mpirun --allow-run-as-root -np 2 ./$(TARGET)

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
	@echo "heFFTe headers:"
	@ls -l /usr/local/include/heffte.h /usr/local/include/heffte_backend_fftw.h
	@echo
	@echo "heFFTe libraries:"
	@ls -l /usr/local/lib/libheffte.so*
	@echo
	@echo "FFTW libraries:"
	@ldconfig -p | grep 'libfftw3'

clean:
	rm -f $(TARGET)
