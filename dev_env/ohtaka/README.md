# Ohtaka / System B Environment

This note records the Ohtaka System B setup used for the oneAPI + Intel MPI
build. It is intended to be reproducible from another Ohtaka account.

The setup below uses:

- Intel oneAPI compiler module
- Intel MPI module
- Intel oneMKL backend for heFFTe
- User-local heFFTe installation under `$HOME/local/heffte-oneapi-mkl`

Do not use `mpirun` on the login node. Use Slurm `srun` inside an allocation or
a batch job.


## Modules

Load the same modules before building heFFTe, building this solver, and running
jobs.

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0
```

Useful checks:

```sh
which icpx
which mpiicpc
which srun
echo $MKLROOT
echo $I_MPI_ROOT
```

Expected paths are similar to:

```text
/opt/intel/oneapi/compiler/2023.0.0/linux/bin/icpx
/opt/intel/oneapi/mpi/2021.8.0/bin/mpiicpc
/usr/bin/srun
/opt/intel/oneapi/mkl/2023.0.0
/opt/intel/oneapi/mpi/2021.8.0
```


## Install heFFTe

Install heFFTe in the account-local prefix. This keeps the setup independent for
each Ohtaka account.

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

mkdir -p $HOME/local/src
mkdir -p $HOME/local/heffte-oneapi-mkl
cd $HOME/local/src

git clone https://github.com/icl-utk-edu/heffte.git
cd heffte
git checkout v2.4.0

mkdir -p build-oneapi-mkl
cd build-oneapi-mkl

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HOME/local/heffte-oneapi-mkl \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_CXX_COMPILER=mpiicpc \
  -DHeffte_ENABLE_MKL=ON \
  -DMKL_ROOT=$MKLROOT \
  -DHeffte_ENABLE_FFTW=OFF \
  -DHeffte_ENABLE_CUDA=OFF \
  -DHeffte_ENABLE_ROCM=OFF \
  -DHeffte_ENABLE_ONEAPI=OFF \
  -DHeffte_ENABLE_AVX=ON \
  -DHeffte_ENABLE_AVX512=OFF

make -j 8
make install
make test_install
```

Set the runtime paths before building or running this solver:

```sh
export HEFFTE_ROOT=$HOME/local/heffte-oneapi-mkl
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

`OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` are recommended for the current flat
MPI runs. `ldd` may show `libmkl_intel_thread.so` and `libiomp5.so` through
heFFTe/MKL, so setting these avoids accidental oversubscription.


## Solver Backend

The Docker development build uses the heFFTe FFTW backend. The Ohtaka build uses
the heFFTe MKL backend. The source should therefore support selecting the backend
with a compile-time macro:

```cpp
#ifdef SHHD_HEFFTE_BACKEND_MKL
     using fft_type = heffte::fft3d_r2c<heffte::backend::mkl>;
#else
     using fft_type = heffte::fft3d_r2c<heffte::backend::fftw>;
#endif
```

This replaces the direct FFTW-only definition in `src/fourier_transform.h`.


## Build The Solver

Use the Ohtaka makefile in the repository root:

```make
CXX := mpiicpc
HEFFTE_ROOT ?= $(HOME)/local/heffte-oneapi-mkl

CXXFLAGS := -std=c++17 -O3 -march=core-avx2 -Wall -Wextra -pedantic \
            -diag-disable=11074,11076 \
            -DSHHD_HEFFTE_BACKEND_MKL \
            -Isrc \
            -I$(HEFFTE_ROOT)/include \
            -I$(MKLROOT)/include

LDFLAGS := -L$(HEFFTE_ROOT)/lib -L$(MKLROOT)/lib/intel64

LDLIBS := -lheffte \
          -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
          -lpthread -lm -ldl
```

Then build:

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-mkl
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8
```

Check linkage:

```sh
ldd src/out.exe | grep -E "heffte|mkl|mpi"
```

Expected libraries include `libheffte.so`, `libmkl_intel_lp64.so`,
`libmkl_sequential.so`, `libmkl_core.so`, and Intel MPI libraries.


## Slurm Run Examples

Interactive test:

```sh
salloc -N 1 -p i8cpu

cd ~/spectral-hohenberg-halperin-dynamics-2d
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-mkl
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun -n 8 -c 1 ./src/out.exe examples/02_incompressible_turbulence/input.script
```

Batch job template:

```sh
#!/bin/sh
#SBATCH -p F1cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -J shhd2d
#SBATCH -o logs/%x-%j.out

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-mkl
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
  ./src/out.exe examples/02_incompressible_turbulence/input.script
```


## Benchmark Conditions

The benchmark below used `examples/02_incompressible_turbulence/input.script`
with:

- `time_evolution srk3/incompressible`
- `fix momentum nonlinear on`
- `fix all noise off`
- `dealias three_halves`
- `measure spec snapshot off`
- `OMP_NUM_THREADS=1`
- `MKL_NUM_THREADS=1`

With `dealias three_halves`, the computational grid is `3/2` times the active
grid in each direction.


## Ohtaka Benchmarks

### Active Grid 128 x 128

`run 10000`, computational grid `192 x 192`:

| MPI ranks | Time [s] | ms/step | Speedup vs 2 ranks | Efficiency vs 2 ranks |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 137.528 | 13.753 | 1.00 | 100% |
| 4 | 78.281 | 7.828 | 1.76 | 87.9% |
| 8 | 41.530 | 4.153 | 3.31 | 82.8% |
| 16 | 24.496 | 2.450 | 5.61 | 70.2% |

Additional short test: `srun -n 1 -c 1`, `run 500` took `11.705 s`
(`23.410 ms/step`). The one-rank number is less reliable as a scaling baseline
because the run length was shorter.

Recommendation for active grid `128 x 128`:

- Use `8` ranks for a good speed/efficiency balance.
- Use `16` ranks when wall-clock time matters more.
- Use `4` ranks when conserving allocation is more important; it is close to
  the best macOS Docker timing measured for this grid.

### Active Grid 512 x 512

`run 100`, computational grid `768 x 768`:

| MPI ranks | Time [s] | ms/step | Speedup vs 1 rank | Efficiency vs 1 rank |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 41.833 | 418.330 | 1.00 | 100% |
| 2 | 27.891 | 278.910 | 1.50 | 75.0% |
| 4 | 13.730 | 137.300 | 3.05 | 76.2% |
| 8 | 8.536 | 85.360 | 4.90 | 61.3% |
| 16 | 6.924 | 69.240 | 6.04 | 37.8% |

`run 500`, computational grid `768 x 768`:

| MPI ranks | Time [s] | ms/step | Speedup vs 4 ranks | Efficiency vs 4 ranks |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 69.637 | 139.274 | 1.00 | 100% |
| 8 | 43.256 | 86.512 | 1.61 | 80.5% |
| 16 | 32.312 | 64.624 | 2.16 | 53.9% |

Recommendation for active grid `512 x 512`:

- Use `8` ranks for a good speed/efficiency balance.
- Use `16` ranks when wall-clock time matters more.
- Use `4` ranks for smaller allocation or queue-friendly test runs.


## macOS Docker Reference

These numbers were measured with the Docker FFTW backend on macOS and are kept
only as a reference. Ohtaka uses oneAPI/MKL backend, so this is not a strict
backend-to-backend comparison.

Active grid `128 x 128`, `run 5000`, snapshot off:

| MPI ranks | Time [s] | ms/step |
| ---: | ---: | ---: |
| 1 | 39.023 | 7.805 |
| 2 | 44.203 | 8.841 |
| 4 | 89.850 | 17.970 |

Active grid `512 x 512`, `run 500`, snapshot off:

| MPI ranks | Time [s] | ms/step |
| ---: | ---: | ---: |
| 1 | 106.230 | 212.460 |
| 2 | 141.199 | 282.398 |
| 4 | 136.714 | 273.428 |

The macOS Docker runs did not show useful MPI scaling. Ohtaka did, even for
`128 x 128`, once the run length was long enough to reduce fixed overhead.


## Notes

- For timing benchmarks, turn snapshot output off. Snapshot I/O should be
  measured separately for production runs.
- Use `srun`, not `mpirun`, on Ohtaka.
- Keep MKL and OpenMP thread counts at `1` for the flat MPI runs above.
- For larger grids such as `1024 x 1024`, benchmark `16`, `32`, and possibly
  more ranks before choosing a production setting.
