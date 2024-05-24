FROM ubuntu:20.04

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        apt-transport-https \
        autoconf \
        automake \
        autotools-dev \
        ca-certificates \
        cmake \
        git \
        libtool \
        m4 \
        make \
        python \
        valgrind \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gfortran-9 && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcov gcov $(which gcov-9) 30 && \
    update-alternatives --install /usr/bin/gfortran gfortran $(which gfortran-9) 30

# LLVM compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        clang-12 \
        clang-format-12 \
        clang-tidy-12 \
        clang-tools-12 \
        libc++-12-dev \
        libc++1-12 \
        libc++abi1-12 \
        libclang-12-dev \
        libclang1-12 \
        liblldb-12-dev \
        libomp-12-dev \
        lld-12 \
        lldb-12 \
        llvm-12 \
        llvm-12-dev \
        llvm-12-runtime && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang $(which clang-12) 30 && \
    update-alternatives --install /usr/bin/clang++ clang++ $(which clang++-12) 30 && \
    update-alternatives --install /usr/bin/clang-format clang-format $(which clang-format-12) 30 && \
    update-alternatives --install /usr/bin/clang-tidy clang-tidy $(which clang-tidy-12) 30 && \
    update-alternatives --install /usr/bin/lldb lldb $(which lldb-12) 30 && \
    update-alternatives --install /usr/bin/llvm-config llvm-config $(which llvm-config-12) 30 && \
    update-alternatives --install /usr/bin/llvm-cov llvm-cov $(which llvm-cov-12) 30

# OFED
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -t focal \
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

# UCX version 1.12.1
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        binutils-dev \
        file \
        libnuma-dev \
        make \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/openucx/ucx/releases/download/v1.12.1/ucx-1.12.1.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/ucx-1.12.1.tar.gz -C /var/tmp -z && \
    cd /var/tmp/ucx-1.12.1 &&  CC=clang CXX=clang++ ./configure --prefix=/usr/local/ucx --disable-assertions --disable-debug --disable-doxygen-doc --disable-logging --disable-params-check --enable-optimizations --without-cuda && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/ucx-1.12.1 /var/tmp/ucx-1.12.1.tar.gz
ENV CPATH=/usr/local/ucx/include:$CPATH \
    LD_LIBRARY_PATH=/usr/local/ucx/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/ucx/lib:$LIBRARY_PATH \
    PATH=/usr/local/ucx/bin:$PATH

# OpenMPI version 4.1.3
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
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://www.open-mpi.org/software/ompi/v4.1/downloads/openmpi-4.1.3.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/openmpi-4.1.3.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/openmpi-4.1.3 &&  CC=clang CXX=clang++ ./configure --prefix=/usr/local/openmpi --disable-mpi-fortran --with-ucx --with-verbs --without-cuda && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/openmpi-4.1.3 /var/tmp/openmpi-4.1.3.tar.bz2
ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/openmpi/bin:$PATH

# HDF5 version 1.12.1
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bzip2 \
        file \
        make \
        wget \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/hdf5-1.12.1.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/hdf5-1.12.1 &&  CC=mpicc CXX=mpicxx F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/hdf5 --disable-fortran --enable-build-mode=production --enable-optimization=high --enable-parallel --with-default-api-version=v110 && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.12.1 /var/tmp/hdf5-1.12.1.tar.bz2
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
    cd /var/tmp/fftw-3.3.10 &&  CC=mpicc CXX=mpicxx F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/fftw --disable-fortran --enable-avx --enable-openmp --enable-sse2 && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/fftw-3.3.10 /var/tmp/fftw-3.3.10.tar.gz
ENV LD_LIBRARY_PATH=/usr/local/fftw/lib:$LD_LIBRARY_PATH

FROM ubuntu:20.04

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        apt-transport-https \
        autoconf \
        automake \
        autotools-dev \
        ca-certificates \
        cmake \
        git \
        libtool \
        m4 \
        make \
        python \
        valgrind \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gfortran-9 && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcov gcov $(which gcov-9) 30 && \
    update-alternatives --install /usr/bin/gfortran gfortran $(which gfortran-9) 30

# LLVM compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        clang-12 \
        clang-format-12 \
        clang-tidy-12 \
        clang-tools-12 \
        libc++-12-dev \
        libc++1-12 \
        libc++abi1-12 \
        libclang-12-dev \
        libclang1-12 \
        liblldb-12-dev \
        libomp-12-dev \
        lld-12 \
        lldb-12 \
        llvm-12 \
        llvm-12-dev \
        llvm-12-runtime && \
    rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang $(which clang-12) 30 && \
    update-alternatives --install /usr/bin/clang++ clang++ $(which clang++-12) 30 && \
    update-alternatives --install /usr/bin/clang-format clang-format $(which clang-format-12) 30 && \
    update-alternatives --install /usr/bin/clang-tidy clang-tidy $(which clang-tidy-12) 30 && \
    update-alternatives --install /usr/bin/lldb lldb $(which lldb-12) 30 && \
    update-alternatives --install /usr/bin/llvm-config llvm-config $(which llvm-config-12) 30 && \
    update-alternatives --install /usr/bin/llvm-cov llvm-cov $(which llvm-cov-12) 30

# GNU compiler runtime
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libgfortran4 \
        libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# LLVM compiler runtime
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libclang1-12 \
        libomp5-12 && \
    rm -rf /var/lib/apt/lists/*

# OFED
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -t focal \
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

# UCX
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libbinutils && \
    rm -rf /var/lib/apt/lists/*
COPY --from=0 /usr/local/ucx /usr/local/ucx
ENV CPATH=/usr/local/ucx/include:$CPATH \
    LD_LIBRARY_PATH=/usr/local/ucx/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/ucx/lib:$LIBRARY_PATH \
    PATH=/usr/local/ucx/bin:$PATH

# OpenMPI
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        hwloc \
        openssh-client && \
    rm -rf /var/lib/apt/lists/*
COPY --from=0 /usr/local/openmpi /usr/local/openmpi
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
ARG FLUPS_DATE=2042-07-19
### h3lpr
RUN --mount=type=ssh git clone git@github.com:van-Rees-Lab/h3lpr.git
RUN cd /h3lpr &&\
    git checkout origin/develop &&\
    echo "CXX = mpic++" > make_arch/make.temp &&\
    echo "CC  = mpicc" >> make_arch/make.temp &&\
    echo "CXXFLAGS := -O3 -std=c++17 -fopenmp -DCOLOR_PROF " >> make_arch/make.temp  &&\
    echo "CCFLAGS := -O3 -std=c99 -fopenmp -DCOLOR_PROF" >> make_arch/make.temp &&\
    echo "LDFLAGS += -fopenmp -lstdc++ -lm" >> make_arch/make.temp &&\
    echo "GTEST_INC = /soft/googletest/include" >> make_arch/make.temp &&\
    echo "GTEST_LIB = /soft/googletest/lib" >> make_arch/make.temp &&\
    ARCH_FILE=make_arch/make.temp PREFIX=/soft/h3lpr make install -j4 &&\
    rm -rf /h3lpr

# add a user
RUN useradd -ms /bin/bash flups
USER flups

# end of file
