# Ohtaka / System B Environment

This note records the Ohtaka System B setup used for the oneAPI + Intel MPI
build. The current recommended configuration uses a user-local FFTW build and
the heFFTe FFTW backend.

The setup below uses:

- Intel oneAPI compiler module
- Intel MPI module
- User-local FFTW installation under `$HOME/local/fftw-oneapi`
- User-local heFFTe installation under `$HOME/local/heffte-oneapi-fftw`

Do not use `mpirun` on the login node. Use Slurm `srun` inside an allocation or
a batch job.


## Modules

Load the same modules before building FFTW, building heFFTe, building this
solver, and running jobs.

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0
```

Useful checks:

```sh
which icpx
which mpiicpc
which srun
echo $I_MPI_ROOT
```

Expected paths are similar to:

```text
/opt/intel/oneapi/compiler/2023.0.0/linux/bin/icpx
/opt/intel/oneapi/mpi/2021.8.0/bin/mpiicpc
/usr/bin/srun
/opt/intel/oneapi/mpi/2021.8.0
```


## Install FFTW

Build both double and float FFTW with Intel MPI. heFFTe links both precisions
when the FFTW backend is enabled.

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export FFTW_ROOT=$HOME/local/fftw-oneapi
export FFTW_VERSION=3.3.10

mkdir -p $HOME/local/src
cd $HOME/local/src

if [ ! -f fftw-$FFTW_VERSION.tar.gz ]; then
  wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
fi

rm -rf fftw-$FFTW_VERSION-double fftw-$FFTW_VERSION-float

tar xzf fftw-$FFTW_VERSION.tar.gz
mv fftw-$FFTW_VERSION fftw-$FFTW_VERSION-double

tar xzf fftw-$FFTW_VERSION.tar.gz
mv fftw-$FFTW_VERSION fftw-$FFTW_VERSION-float

export CC=mpiicc
export MPICC=mpiicc

cd $HOME/local/src/fftw-$FFTW_VERSION-double

./configure \
  --prefix=$FFTW_ROOT \
  --enable-shared \
  --enable-mpi \
  --enable-sse2 \
  --enable-avx \
  --enable-avx2 \
  --disable-fortran \
  CFLAGS="-O3 -march=core-avx2"
make -j 8
make install

cd $HOME/local/src/fftw-$FFTW_VERSION-float

./configure \
  --prefix=$FFTW_ROOT \
  --enable-shared \
  --enable-mpi \
  --enable-float \
  --enable-sse2 \
  --enable-avx \
  --enable-avx2 \
  --disable-fortran \
  CFLAGS="-O3 -march=core-avx2"
make -j 8
make install
```

Check that the required libraries exist:

```sh
ls -lh \
  $FFTW_ROOT/lib/libfftw3.so \
  $FFTW_ROOT/lib/libfftw3_mpi.so \
  $FFTW_ROOT/lib/libfftw3f.so \
  $FFTW_ROOT/lib/libfftw3f_mpi.so
```


## Install heFFTe

Install heFFTe in an account-local prefix and use the FFTW backend.

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export FFTW_ROOT=$HOME/local/fftw-oneapi
export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export LD_LIBRARY_PATH=$FFTW_ROOT/lib:$LD_LIBRARY_PATH

mkdir -p $HOME/local/src
cd $HOME/local/src

git clone --depth 1 --branch v2.4.0 https://github.com/icl-utk-edu/heffte.git heffte-fftw
cd heffte-fftw

mkdir -p build-oneapi-fftw
cd build-oneapi-fftw

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HEFFTE_ROOT \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_CXX_COMPILER=mpiicpc \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=core-avx2" \
  -DCMAKE_PREFIX_PATH=$FFTW_ROOT \
  -DFFTW_ROOT=$FFTW_ROOT \
  -DHeffte_ENABLE_FFTW=ON \
  -DHeffte_ENABLE_MKL=OFF \
  -DHeffte_ENABLE_CUDA=OFF \
  -DHeffte_ENABLE_ROCM=OFF \
  -DHeffte_ENABLE_ONEAPI=OFF \
  -DHeffte_ENABLE_AVX=ON \
  -DHeffte_ENABLE_AVX512=OFF \
  -DHeffte_ENABLE_TESTING=OFF

cmake --build . --target Heffte --parallel 8
cmake --install .
```

Check linkage:

```sh
ldd $HEFFTE_ROOT/lib/libheffte.so | grep -E "fftw|mkl|mpi"
```

Expected: FFTW and Intel MPI libraries appear. MKL libraries should not appear
as FFT backend dependencies.


## Build The Solver

Use the Ohtaka makefile in the repository root:

```sh
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE

make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8
```

Check linkage:

```sh
ldd src/out.exe | grep -E "heffte|fftw|mkl|mpi"
```

Expected: `libheffte.so`, `libfftw3_mpi.so`, `libfftw3.so`,
`libfftw3f_mpi.so`, `libfftw3f.so`, and Intel MPI libraries. MKL libraries
should not appear.


## Slurm Run Examples

Interactive test:

```sh
salloc -N 1 -p i8cpu

cd ~/spectral-hohenberg-halperin-dynamics-2d
module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE

srun --cpu-bind=cores -n 8 -c 1 \
  ./src/out.exe examples/02_incompressible_turbulence/input.script
```

Batch job template for reserving one full `i8cpu` node but running a selected
number of MPI ranks:

```sh
#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 00:30:00
#SBATCH -J shhd2d
#SBATCH -o shhd2d-%j.out
#SBATCH -e shhd2d-%j.err

set -eu

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE

export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core

# Reserve the full node, but run 16 MPI ranks.
srun --cpu-bind=cores --distribution=block:block \
  -n 16 -c 1 \
  ./src/out.exe examples/02_incompressible_turbulence/input.script
```

For the current x-slab R2C decomposition, the number of MPI ranks must satisfy
`size <= Nx/2 + 1`. For `Nx = 128`, ranks up to `64` are safe; `128` ranks are
not valid.


### Yokota Green-Kubo Linear Samples

`examples/05_yokota_green_kubo` has a one-node sample job for the linear
Yokota Green-Kubo check. Generate the per-sample input scripts before
submitting the Slurm job:

```sh
cd ~/spectral-hohenberg-halperin-dynamics-2d

python3 examples/05_yokota_green_kubo/prepare_ohtaka_linear_inputs.py \
  --run-id prepared \
  --samples 128

sbatch examples/05_yokota_green_kubo/job_ohtaka_linear_samples.sh
```

The job reserves one `i8cpu` node and launches 128 independent one-rank runs
with one `srun` per sample. It intentionally does not use a Slurm array, so it
avoids array submission limits and keeps all samples inside one batch job.

The important launch options are:

```sh
srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n 1 -c 1 -N 1 ...
```

`--cpu-bind=cores` is important on Ohtaka for avoiding the large run-to-run
speed variation seen when many one-rank samples are launched concurrently.

To use another run id, generate matching inputs and export the same id when
submitting:

```sh
python3 examples/05_yokota_green_kubo/prepare_ohtaka_linear_inputs.py \
  --run-id linear_test \
  --samples 128

YKGK_RUN_ID=linear_test \
  sbatch examples/05_yokota_green_kubo/job_ohtaka_linear_samples.sh
```

The generated inputs are written under
`examples/05_yokota_green_kubo/runs/ohtaka_linear_<run-id>/`, and outputs are
written under
`examples/05_yokota_green_kubo/results/ohtaka_linear_<run-id>/`.


## Benchmark Conditions

The benchmark below used `examples/02_incompressible_turbulence/input.script`
with:

- active grid `128 x 128`
- `run 5000`
- one full `i8cpu` node reserved
- flat MPI runs with `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1`
- `srun --cpu-bind=cores --distribution=block:block -n <ranks> -c 1`


## Ohtaka Benchmarks

### Active Grid 128 x 128

| MPI ranks | Time [s] | Speedup vs 1 rank | Efficiency vs 1 rank |
| ---: | ---: | ---: | ---: |
| 1 | 53.227 | 1.00 | 100% |
| 2 | 31.393 | 1.70 | 84.8% |
| 4 | 17.254 | 3.09 | 77.1% |
| 8 | 10.061 | 5.29 | 66.1% |
| 16 | 5.976 | 8.90 | 55.7% |

Recommendation for active grid `128 x 128`:

- Use `8` ranks for a good speed/efficiency balance.
- Use `16` ranks when wall-clock time matters more.
- There is little need to go beyond `16` ranks for this grid unless the input
  physics changes the cost balance.


## macOS Docker Reference

These numbers were measured with the Docker FFTW backend on macOS and are kept
only as a reference.

Active grid `128 x 128`, `run 5000`, snapshot off:

| MPI ranks | Time [s] | ms/step |
| ---: | ---: | ---: |
| 1 | 39.023 | 7.805 |
| 2 | 44.203 | 8.841 |
| 4 | 89.850 | 17.970 |

The macOS Docker runs did not show useful MPI scaling. Ohtaka did, once the
oneAPI + Intel MPI + FFTW backend environment was used.


## Cleanup

After installation, the build directories under `$HOME/local/src` can be
removed if disk usage matters:

```sh
rm -rf $HOME/local/src
```

Keep these installed prefixes:

```text
$HOME/local/fftw-oneapi
$HOME/local/heffte-oneapi-fftw
```


## Notes

- Use `srun`, not `mpirun`, on Ohtaka.
- Keep OpenMP and MKL thread counts at `1` for the flat MPI runs above.
- For timing benchmarks, turn snapshot output off. Snapshot I/O should be
  measured separately for production runs.
- For larger grids such as `512 x 512` and `1024 x 1024`, benchmark `8`, `16`,
  `32`, and possibly more ranks before choosing a production setting.
