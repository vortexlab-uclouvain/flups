FROM gilponcelet/gpu-benchmarks


USER root
WORKDIR /

# RUN echo "nameserver 8.8.4.4" >> /etc/resolv.conf
#RUN apt update
#RUN apt upgrade -y
RUN apt -y update && apt install -y \
make git g++ libhdf5-mpi-dev linux-libc-dev vim wget unzip cmake
#valgrind

#METIS:
RUN wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
RUN tar -xvf metis-5.1.0.tar.gz
RUN rm metis-5.1.0.tar.gz
RUN cd /metis-5.1.0 && make config cc=gcc && make && make install

#Google test
RUN git clone https://github.com/google/googletest.git
RUN cd /googletest &&\
    cmake . -DCMAKE_INSTALL_PREFIX=/soft/googletest &&\
    make install -j
RUN rm -rf /googletest

# h3lpr
COPY h3lpr /h3lpr
RUN cd /h3lpr &&\
    echo "CXX = mpic++" > make_arch/make.temp &&\
    echo "CC = mpicc" >> make_arch/make.temp &&\
    echo "CXXFLAGS := -O3 -std=c++17 -fopenmp -DCOLOR_PROF " >> make_arch/make.temp  &&\
    echo "CCFLAGS := -O3 -std=c99 -fopenmp -DCOLOR_PROF" >> make_arch/make.temp &&\
    echo "LDFLAGS += -fopenmp -lstdc++ -lm" >> make_arch/make.temp &&\
    echo "GTEST_INC = /soft/googletest/include" >> make_arch/make.temp &&\
    echo "GTEST_LIB = /soft/googletest/lib" >> make_arch/make.temp &&\
    ARCH_FILE=make_arch/make.temp PREFIX=/soft/h3lpr make install -j4 &&\
    rm -rf /h3lpr


# VOLUME ["/etc/gitlab-runner"]

# RUN useradd -ms /bin/bash validator
USER dev
WORKDIR /home/dev
