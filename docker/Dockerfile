FROM ubuntu:18.04

RUN apt update
RUN apt install -y \
make git g++ libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin python3.7

VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
WORKDIR /home/validator
