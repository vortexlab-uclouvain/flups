FROM ubuntu:18.04


RUN echo "nameserver 8.8.4.4" >> /etc/resolv.conf
RUN apt update
#RUN apt upgrade -y

RUN apt install -y \
<<<<<<< HEAD
make git g++ libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin python3.7 linux-libc-dev valgrind
=======
make git g++ vim libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin python3.7
>>>>>>> dev_api_c

VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
WORKDIR /home/validator
