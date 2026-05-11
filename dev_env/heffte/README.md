# heFFTe Environment Test

This directory contains a small direct test of heFFTe with the FFTW backend.

The test uses a 4 by 4 by 4 complex grid split across two MPI ranks in the
third direction. It performs:

1. A forward 3D FFT with heFFTe.
2. A backward 3D FFT with heFFTe.
3. A comparison against the original local input after dividing by the total
   number of grid points.

The test is intentionally small and uses a traditional Makefile so that the
include paths and link libraries stay visible.

## Manual Compile

From inside the Docker container, at the repository root:

```bash
mpic++ -std=c++14 \
  -I/usr/local/include \
  dev_env/heffte/check_heffte_fftw.cpp \
  -L/usr/local/lib \
  -lheffte \
  -lfftw3_mpi -lfftw3 \
  -lfftw3f_mpi -lfftw3f \
  -lfftw3l_mpi -lfftw3l \
  -o dev_env/heffte/out.exe
```

Then run:

```bash
mpirun --allow-run-as-root -np 2 dev_env/heffte/out.exe
```

Expected output:

```text
heFFTe FFTW check OK with 2 ranks
```

## Makefile

From inside the Docker container, at the repository root:

```bash
make -C dev_env/heffte
make -C dev_env/heffte run
```

To inspect the include and library paths:

```bash
make -C dev_env/heffte paths
```

To clean:

```bash
make -C dev_env/heffte clean
```

## Notes

The extra FFTW libraries are linked because this heFFTe build exposes support
for double, float, and long-double FFTW backends:

```text
-lfftw3_mpi -lfftw3
-lfftw3f_mpi -lfftw3f
-lfftw3l_mpi -lfftw3l
```

The heFFTe headers and library are installed in the Docker image under:

```text
/usr/local/include/heffte.h
/usr/local/lib/libheffte.so
```
