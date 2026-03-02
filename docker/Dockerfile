FROM ubuntu:24.04

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        apt-transport-https \
        autoconf \
        automake \
        autotools-dev \
        ca-certificates \
        cmake \
        git \
        libstdc++-14-dev \
        libtool \
        libucx-dev \
        m4 \
        make \
        valgrind \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gfortran-14 && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcov gcov $(which gcov-14) 30 && \
    update-alternatives --install /usr/bin/gfortran gfortran $(which gfortran-14) 30

# LLVM compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        clang-18 \
        clang-format-18 \
        clang-tidy-18 \
        clang-tools-18 \
        libc++-18-dev \
        libc++1-18 \
        libc++abi1-18 \
        libclang-18-dev \
        libclang1-18 \
        liblldb-18-dev \
        libomp-18-dev \
        lld-18 \
        lldb-18 \
        llvm-18 \
        llvm-18-dev \
        llvm-18-runtime && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang $(which clang-18) 30 && \
    update-alternatives --install /usr/bin/clang++ clang++ $(which clang++-18) 30 && \
    update-alternatives --install /usr/bin/clang-format clang-format $(which clang-format-18) 30 && \
    update-alternatives --install /usr/bin/clang-tidy clang-tidy $(which clang-tidy-18) 30 && \
    update-alternatives --install /usr/bin/lldb lldb $(which lldb-18) 30 && \
    update-alternatives --install /usr/bin/llvm-config llvm-config $(which llvm-config-18) 30 && \
    update-alternatives --install /usr/bin/llvm-cov llvm-cov $(which llvm-cov-18) 30

# OFED
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -t noble \
        dapl2-utils \
        ibutils \
        ibverbs-providers \
        ibverbs-utils \
        infiniband-diags \
        libdapl-dev \
        libdapl2 \
        libibmad-dev \
        libibmad5 \
        libibverbs-dev \
        libibverbs1 \
        librdmacm-dev \
        librdmacm1 \
        rdmacm-utils && \
    rm -rf /var/lib/apt/lists/*

# OpenMPI version 5.0.9
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bzip2 \
        file \
        hwloc \
        libnuma-dev \
        make \
        openssh-client \
        perl \
        tar \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://www.open-mpi.org/software/ompi/v5.0/downloads/openmpi-5.0.9.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/openmpi-5.0.9.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/openmpi-5.0.9 &&  CC=clang CXX=clang++ ./configure --prefix=/usr/local/openmpi --disable-mpi-fortran --with-ucx --with-verbs --without-cuda && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/openmpi-5.0.9 /var/tmp/openmpi-5.0.9.tar.bz2
ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/openmpi/bin:$PATH

# HDF5 version 1.14.6
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bzip2 \
        file \
        make \
        wget \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/hdf5-1.14.6.tar.gz -C /var/tmp -z && \
    cd /var/tmp/hdf5-1.14.6 &&  CC=mpicc CXX=mpicxx F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/hdf5 --disable-fortran --enable-build-mode=production --enable-optimization=high --enable-parallel --with-default-api-version=v110 && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.14.6 /var/tmp/hdf5-1.14.6.tar.gz
ENV CPATH=/usr/local/hdf5/include:$CPATH \
    HDF5_DIR=/usr/local/hdf5 \
    LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/hdf5/lib:$LIBRARY_PATH \
    PATH=/usr/local/hdf5/bin:$PATH

# FFTW version 3.3.10
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        file \
        make \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp ftp://ftp.fftw.org/pub/fftw/fftw-3.3.10.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/fftw-3.3.10.tar.gz -C /var/tmp -z && \
    cd /var/tmp/fftw-3.3.10 &&  CC=mpicc CXX=mpicxx F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/fftw --disable-fortran --enable-openmp --enable-shared --enable-static && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/fftw-3.3.10 /var/tmp/fftw-3.3.10.tar.gz
ENV LD_LIBRARY_PATH=/usr/local/fftw/lib:$LD_LIBRARY_PATH

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python && \
    pip3 install numpy --break-system-packages && \
    MPICC=/usr/local/openmpi/bin/mpicc pip3 install mpi4py --break-system-packages --no-cache-dir
FROM ubuntu:24.04

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        apt-transport-https \
        autoconf \
        automake \
        autotools-dev \
        ca-certificates \
        cmake \
        git \
        libstdc++-14-dev \
        libtool \
        libucx-dev \
        m4 \
        make \
        valgrind \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*

# LLVM compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        clang-18 \
        clang-format-18 \
        clang-tidy-18 \
        clang-tools-18 \
        libc++-18-dev \
        libc++1-18 \
        libc++abi1-18 \
        libclang-18-dev \
        libclang1-18 \
        liblldb-18-dev \
        libomp-18-dev \
        lld-18 \
        lldb-18 \
        llvm-18 \
        llvm-18-dev \
        llvm-18-runtime && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang $(which clang-18) 30 && \
    update-alternatives --install /usr/bin/clang++ clang++ $(which clang++-18) 30 && \
    update-alternatives --install /usr/bin/clang-format clang-format $(which clang-format-18) 30 && \
    update-alternatives --install /usr/bin/clang-tidy clang-tidy $(which clang-tidy-18) 30 && \
    update-alternatives --install /usr/bin/lldb lldb $(which lldb-18) 30 && \
    update-alternatives --install /usr/bin/llvm-config llvm-config $(which llvm-config-18) 30 && \
    update-alternatives --install /usr/bin/llvm-cov llvm-cov $(which llvm-cov-18) 30

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gfortran-14 && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcov gcov $(which gcov-14) 30 && \
    update-alternatives --install /usr/bin/gfortran gfortran $(which gfortran-14) 30

# OpenMPI version 5.0.9
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bzip2 \
        file \
        hwloc \
        libnuma-dev \
        make \
        openssh-client \
        perl \
        tar \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://www.open-mpi.org/software/ompi/v5.0/downloads/openmpi-5.0.9.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/openmpi-5.0.9.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/openmpi-5.0.9 &&  CC=clang CXX=clang++ ./configure --prefix=/usr/local/openmpi --disable-mpi-fortran --with-ucx --with-verbs --without-cuda && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/openmpi-5.0.9 /var/tmp/openmpi-5.0.9.tar.bz2
ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/openmpi/bin:$PATH

# HDF5
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        zlib1g && \
    rm -rf /var/lib/apt/lists/*
COPY --from=0 /usr/local/hdf5 /usr/local/hdf5
ENV CPATH=/usr/local/hdf5/include:$CPATH \
    HDF5_DIR=/usr/local/hdf5 \
    LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/hdf5/lib:$LIBRARY_PATH \
    PATH=/usr/local/hdf5/bin:$PATH

# FFTW
COPY --from=0 /usr/local/fftw /usr/local/fftw
ENV LD_LIBRARY_PATH=/usr/local/fftw/lib:$LD_LIBRARY_PATH

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3 \
        python3-dev \
        python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN rm /usr/lib/python3.*/EXTERNALLY-MANAGED && \
    ln -s /usr/bin/python3 /usr/bin/python || true && \
    pip3 install numpy --break-system-packages && \
    MPICC=/usr/local/openmpi/bin/mpicc pip3 install mpi4py --break-system-packages --force-reinstall --no-cache-dir



## this file is added at the end of the dockerfile

#Google test
RUN git clone https://github.com/google/googletest.git
RUN cd /googletest &&\
    cmake . -DCMAKE_INSTALL_PREFIX=/soft/googletest &&\
    make install -j
RUN rm -rf /googletest

#-------------------------------------------------------------------------------
RUN apt-get -y update && apt-get install -y openssh-client git
RUN mkdir -p -m 0600 ~/.ssh && ssh-keyscan git.immc.ucl.ac.be github.com >> ~/.ssh/known_hosts

#-------------------------------------------------------------------------------
ARG FLUPS_DATE=2026-02-03

### h3lpr
RUN --mount=type=ssh git clone git@github.com:van-Rees-Lab/h3lpr.git
RUN cd /h3lpr &&\
    git checkout origin/develop &&\
    echo "CXX = mpic++" > make_arch/make.temp &&\
    echo "CC  = mpicc" >> make_arch/make.temp &&\
    echo "CXXFLAGS := -O3 -std=c++17 -fopenmp -DCOLOR_PROF" >> make_arch/make.temp  &&\
    echo "CCFLAGS := -O3 -std=c99 -fopenmp -DCOLOR_PROF" >> make_arch/make.temp &&\
    echo "LDFLAGS += -fopenmp -lm" >> make_arch/make.temp &&\
    echo "GTEST_INC = /soft/googletest/include" >> make_arch/make.temp &&\
    echo "GTEST_LIB = /soft/googletest/lib" >> make_arch/make.temp &&\
    ARCH_FILE=make_arch/make.temp PREFIX=/soft/h3lpr make install -j4 &&\
    rm -rf /h3lpr

# add a user
RUN useradd -ms /bin/bash flups 
USER flups
RUN echo "export ARCH_FILE=make_arch/make.docker_clang" >> ~/.bashrc

# end of file
