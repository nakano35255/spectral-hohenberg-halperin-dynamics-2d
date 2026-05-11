FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    git \
    cmake \
    openmpi-bin \
    libopenmpi-dev \
    libfftw3-dev \
    libfftw3-mpi-dev

RUN git clone --depth 1 --branch v2.4.1 https://github.com/icl-utk-edu/heffte.git /tmp/heffte \
    && cmake -S /tmp/heffte -B /tmp/heffte/build \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=/usr/local \
        -D BUILD_SHARED_LIBS=ON \
        -D Heffte_ENABLE_FFTW=ON \
        -D MPI_CXX_COMPILER=mpic++ \
    && cmake --build /tmp/heffte/build -j4 \
    && cmake --install /tmp/heffte/build \
    && ldconfig \
    && rm -rf /tmp/heffte