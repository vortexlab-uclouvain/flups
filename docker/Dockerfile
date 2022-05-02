FROM ubuntu:18.04


# RUN echo "nameserver 8.8.4.4" >> /etc/resolv.conf
#RUN apt update
#RUN apt upgrade -y
RUN apt -y update && apt install -y \
make git g++ libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin python3.7 linux-libc-dev vim wget unzip cmake
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


VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
WORKDIR /home/validator
