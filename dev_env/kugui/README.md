# Kugui / System C Environment

This note records the Kugui System C setup used for the Intel compiler + Intel
MPI build. The recommended configuration mirrors the Ohtaka setup: build FFTW
under the user account and build heFFTe with the FFTW backend.

The setup below uses:

- Intel compiler module
- Intel MPI module
- User-local FFTW installation under `$HOME/local/fftw-oneapi`
- User-local heFFTe installation under `$HOME/local/heffte-oneapi-fftw`
- PBS jobs for execution

Do not use `mpirun` for production runs on the login node. Use PBS `qsub` and
run `mpirun` inside the batch job.


## Login

```sh
ssh i001900@kugui.issp.u-tokyo.ac.jp
```

The login node used during setup was `kugui1`.


## Modules

Load the same modules before building FFTW, building heFFTe, building this
solver, and running jobs.

```sh
module purge
module load intel/2022.2.1 intel-mpi/2021.7.1
```

Useful checks:

```sh
which icpx
which mpiicpc
which mpirun
echo $I_MPI_ROOT
echo $MKLROOT
```

Expected paths are similar to:

```text
/home/app/oneapi/compiler/2022.2.1/linux/bin/icpx
/home/app/oneapi/mpi/2021.7.1/bin/mpiicpc
/home/app/oneapi/mpi/2021.7.1/bin/mpirun
/home/app/oneapi/mpi/2021.7.1
/home/app/oneapi/mkl/2022.2.1
```

`mpiicc` and `mpiicpc` use Intel classic compilers on this system, so builds may
print Intel compiler deprecation remarks. Those remarks are expected.


## Install FFTW

Build both double and float FFTW with Intel MPI. heFFTe links both precisions
when the FFTW backend is enabled.

```sh
module purge
module load intel/2022.2.1 intel-mpi/2021.7.1

export FFTW_ROOT=$HOME/local/fftw-oneapi
export FFTW_VERSION=3.3.10
export SRC_ROOT=$HOME/local/src

mkdir -p "$SRC_ROOT" "$FFTW_ROOT"
cd "$SRC_ROOT"

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

cd "$SRC_ROOT/fftw-$FFTW_VERSION-double"

./configure \
  --prefix="$FFTW_ROOT" \
  --enable-shared \
  --enable-mpi \
  --enable-sse2 \
  --enable-avx \
  --enable-avx2 \
  --disable-fortran \
  CFLAGS="-O3 -march=core-avx2"
make -j 8
make install

cd "$SRC_ROOT/fftw-$FFTW_VERSION-float"

./configure \
  --prefix="$FFTW_ROOT" \
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
  "$FFTW_ROOT/lib/libfftw3.so" \
  "$FFTW_ROOT/lib/libfftw3_mpi.so" \
  "$FFTW_ROOT/lib/libfftw3f.so" \
  "$FFTW_ROOT/lib/libfftw3f_mpi.so"
```


## Install heFFTe

Install heFFTe in an account-local prefix and use the FFTW backend. MKL is not
used as the FFT backend in this setup.

```sh
module purge
module load intel/2022.2.1 intel-mpi/2021.7.1

export FFTW_ROOT=$HOME/local/fftw-oneapi
export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export LD_LIBRARY_PATH=$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}

mkdir -p $HOME/local/src
cd $HOME/local/src

git clone --depth 1 --branch v2.4.0 https://github.com/icl-utk-edu/heffte.git heffte-fftw
cd heffte-fftw

rm -rf build-oneapi-fftw
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
ldd "$HEFFTE_ROOT/lib/libheffte.so" | grep -E "fftw|mkl|mpi"
```

Expected: FFTW and Intel MPI libraries appear. MKL libraries should not appear
as FFT backend dependencies.


## Clone The Repository

The Kugui account used during setup did not have a GitHub SSH key configured.
HTTPS clone worked:

```sh
cd $HOME
git clone https://github.com/nakano35255/spectral-hohenberg-halperin-dynamics-2d.git
cd spectral-hohenberg-halperin-dynamics-2d
```

If GitHub SSH access is configured later, the SSH remote is also fine:

```sh
git clone git@github.com:nakano35255/spectral-hohenberg-halperin-dynamics-2d.git
```


## Build The Solver

The current Kugui setup can use the same oneAPI/Intel-MPI makefile as Ohtaka:

```sh
module purge
module load intel/2022.2.1 intel-mpi/2021.7.1

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE

make -f Makefile.ohtaka yes-PASSIVE-SCALAR
make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8 all
```

Check linkage:

```sh
ldd src/out.exe | grep -E "heffte|fftw|mkl|mpi"
```

Expected: `libheffte.so`, `libfftw3_mpi.so`, `libfftw3.so`,
`libfftw3f_mpi.so`, `libfftw3f.so`, and Intel MPI libraries. MKL libraries
should not appear.


## Optional Environment File

For repeated use, create a small environment file:

```sh
cat > "$HOME/local/env_shhd_kugui.sh" <<'ENV'
module purge
module load intel/2022.2.1 intel-mpi/2021.7.1
export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE
ENV
```

Then use:

```sh
source "$HOME/local/env_shhd_kugui.sh"
```


## PBS Queues

Kugui uses PBS. Queue limits observed during setup:

```text
i2cpu   1-2 nodes, max 256 cores, 00:30:00
B1cpu   1 node,    max 128 cores, 12:00:00
F1cpu   1 node,    max 128 cores, 24:00:00
L1cpu   1 node,    max 128 cores, 120:00:00
L4cpu   2-4 nodes, max 512 cores, 120:00:00
L16cpu  5-16 nodes, max 2048 cores, 120:00:00
```

Check current queue status with:

```sh
qstat -Q
qstat -Qf L1cpu
qstat -Qf L4cpu
```


## PBS Run Example

A minimal PBS job for a one-node run:

```sh
#!/bin/bash
#PBS -q L1cpu
#PBS -N shhd_run
#PBS -l select=1:ncpus=128:mpiprocs=128:ompthreads=1
#PBS -l walltime=120:00:00
#PBS -j oe

set -euo pipefail

source "$HOME/local/env_shhd_kugui.sh"
cd "$HOME/spectral-hohenberg-halperin-dynamics-2d"

np=$(wc -l < "$PBS_NODEFILE")
echo "PBS_JOBID=${PBS_JOBID:-}"
echo "PBS_NODEFILE=$PBS_NODEFILE"
echo "np=$np"

mpirun -np "$np" -machinefile "$PBS_NODEFILE" \
  ./src/out.exe examples/04_ness_uniform_gradient/cascade/input.script
```

For a four-node run, use `L4cpu` and a matching `select` line:

```sh
#PBS -q L4cpu
#PBS -l select=4:ncpus=128:mpiprocs=128:ompthreads=1
#PBS -l walltime=120:00:00
```

Submit and inspect:

```sh
qsub job_kugui.pbs
qstat -u "$USER"
```


## Smoke Test

The setup was verified on 2026-06-18 with a 2-rank PBS job on `i2cpu`.

```text
PBS_JOBID=824339.kugui-pbs
np=2
Run 1 finished: global step 2 time 2
Simulation Finished. Total time: 0.003 s
```

The tested executable linked to:

```text
libheffte.so
libfftw3_mpi.so
libfftw3.so
libfftw3f_mpi.so
libfftw3f.so
libmpi.so
```


## Notes

- This setup intentionally uses the heFFTe FFTW backend, not the MKL backend.
- `mpiicc` and `mpiicpc` emit Intel classic compiler deprecation remarks on
  Kugui; those messages are not build failures.
- When using the dealiased `three_halves` grid, avoid requesting more MPI ranks
  than the available spectral slab count for the chosen grid.
