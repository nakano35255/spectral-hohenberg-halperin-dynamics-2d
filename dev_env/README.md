# Development Environment

This directory keeps small, explicit tests for the Docker-based development
environment.

## Directories

- `docker/`: Docker, VS Code Dev Containers, MPI, and FFTW setup notes.
- `heffte/`: A direct heFFTe + FFTW + MPI test using a traditional Makefile.

## Quick Checks

From inside the Docker container, at the repository root:

```bash
make -C dev_env/docker clean run
make -C dev_env/heffte clean run
```

Expected output includes:

```text
FFTW MPI check OK with 4 ranks
heFFTe FFTW check OK with 2 ranks
```

For Docker and VS Code setup details, see `dev_env/docker/README.md`.
For the heFFTe link command and test details, see `dev_env/heffte/README.md`.
