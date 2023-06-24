FROM ubuntu:22.04

RUN apt-get update

RUN apt install -y tzdata
ENV TZ=Asia/Tokyo

RUN apt install -y wget \
  gpg \
  g++ \
  git \
  clang-format

# cmake 3.26
RUN cd $HOME && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null && \
    echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list > /dev/null && \
    apt-get update && \
    apt install -y cmake

# boost 1.82.0
RUN cd $HOME && \
    git clone https://github.com/boostorg/boost.git && \
    cd boost && \
    git checkout boost-1.82.0 && \
    git submodule update --init && \
    ./bootstrap.sh && \
    ./b2 toolset=gcc-11 --prefix=/usr/local -j9 && \
    ./b2 install toolset=gcc-11 --prefix=/usr/local -j9
