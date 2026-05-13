# Docker Environment

This directory records how to set up and test the Docker-based development
environment for `spectral-fluctuating-isothermal-fluid`.

The goal is to keep the host machine clean while using the same Linux-based
C++/MPI/FFTW/heFFTe environment on different Macs.

## Requirements on the Mac

Install these on the host machine:

- Docker Desktop
- Visual Studio Code
- VS Code extension: Dev Containers

Check that Docker is available:

```bash
docker --version
docker run hello-world
```

If `Hello from Docker!` appears, Docker is working.

## Build the Docker Image

From the repository root:

```bash
cd /Users/nakano/github/spectral-fluctuating-isothermal-fluid
docker build -t spectral-fluctuating-isothermal-fluid .
```

This reads the repository `Dockerfile` and builds an image containing:

- Ubuntu 24.04
- C++ compiler and build tools
- OpenMPI
- FFTW3
- FFTW3 MPI support
- heFFTe with the FFTW backend

Check the image:

```bash
docker images
```

## Run a Shell in the Container

From the repository root:

```bash
docker run --rm -it \
  -v "$PWD":/work \
  -w /work \
  spectral-fluctuating-isothermal-fluid \
  bash
```

Meaning:

- `--rm`: remove the container after it exits
- `-it`: open an interactive terminal
- `-v "$PWD":/work`: mount the current Mac directory at `/work` in the container
- `-w /work`: start in `/work`

Files created under `/work` are actually files in the Mac repository, so they
remain after the container exits.

Files created elsewhere in the container, such as `/tmp` or `/root`, disappear
when the container exits because of `--rm`.

## Check Installed Tools

Inside the container:

```bash
g++ --version
mpic++ --version
mpirun --version
ls /usr/include/fftw3.h
ls /usr/include/fftw3-mpi.h
ls /usr/local/include/heffte.h
ls /usr/local/lib/libheffte.so
```

Check where MPI, FFTW, and heFFTe are installed:

```bash
which g++
which mpic++
which mpirun
ls -l /usr/include/fftw3.h /usr/include/fftw3-mpi.h
ls -l /usr/local/include/heffte.h
ls -l /usr/local/lib/libheffte.so
ldconfig -p | grep 'libfftw3'
mpic++ --showme:compile
mpic++ --showme:link
```

On an Apple Silicon Mac, the Docker container is usually an ARM64 Linux
environment, so FFTW libraries are typically under:

```text
/lib/aarch64-linux-gnu/
```

These are paths inside the Docker container, not paths on the Mac host.

## Run the MPI + FFTW Test by Hand

For learning, first compile the test program directly.

From inside the container, at the repository root:

```bash
mpic++ -std=c++14 \
  -I/usr/include \
  dev_env/docker/check_mpi_fftw.cpp \
  -L/lib/aarch64-linux-gnu \
  -lfftw3_mpi \
  -lfftw3 \
  -o dev_env/docker/out.exe
```

Then run it with 4 MPI processes:

```bash
mpirun --allow-run-as-root -np 4 dev_env/docker/out.exe
```

Expected output includes:

```text
FFTW MPI check OK with 4 ranks
```

What the compile command means:

- `mpic++`: MPI-aware C++ compiler wrapper
- `-std=c++14`: use the C++14 language standard
- `-I/usr/include`: search this directory for headers such as `fftw3.h`
- `dev_env/docker/check_mpi_fftw.cpp`: source file to compile
- `-L/lib/aarch64-linux-gnu`: search this directory for libraries
- `-lfftw3_mpi`: link `libfftw3_mpi.so`
- `-lfftw3`: link `libfftw3.so`
- `-o dev_env/docker/out.exe`: write the executable here

In this Docker image, `/usr/include` and `/lib/aarch64-linux-gnu` are standard
compiler/linker search paths, so the shorter command also works:

```bash
mpic++ -std=c++14 dev_env/docker/check_mpi_fftw.cpp \
  -lfftw3_mpi -lfftw3 \
  -o dev_env/docker/out.exe
```

The longer command is useful for learning because it makes the include and
library paths explicit.

## Run the MPI + FFTW Test with Make

From inside the container, at the repository root:

```bash
make -C dev_env/docker
make -C dev_env/docker run
```

To inspect the paths used by the Makefile:

```bash
make -C dev_env/docker paths
```

To clean the test executable:

```bash
make -C dev_env/docker clean
```

The `--allow-run-as-root` option is needed because the default Docker container
user is root. This is normal for this simple development image.

## Using VS Code Dev Containers

Open the repository in VS Code:

```bash
cd /Users/nakano/github/spectral-fluctuating-isothermal-fluid
code .
```

Then in VS Code:

```text
Cmd + Shift + P
Dev Containers: Reopen in Container
```

VS Code will build or reuse the Docker image, start a container, and mount the
repository inside it. The VS Code window stays on the Mac, but the terminal,
compiler, MPI, FFTW, and heFFTe commands run inside the Linux container.

After reopening in the container, the VS Code terminal should be inside the
container. Check with:

```bash
pwd
g++ --version
mpic++ --version
ls /usr/local/include/heffte.h
```

If the Dockerfile or `.devcontainer/devcontainer.json` is changed, use:

```text
Cmd + Shift + P
Dev Containers: Rebuild and Reopen in Container
```

For normal daily work, `Dev Containers: Reopen in Container` is enough.

## Notes for Another Mac

On a different Mac:

1. Install Docker Desktop.
2. Install VS Code.
3. Install the VS Code Dev Containers extension.
4. Clone this repository.
5. Run `docker build -t spectral-fluctuating-isothermal-fluid .`.
6. Open the repository in VS Code.
7. Use `Dev Containers: Reopen in Container`.
8. Run the environment tests.

This should reproduce the same C++/MPI/FFTW/heFFTe development environment.
