cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_COMPILER /system/apps/ubuntu/20.04-202304/nvhpc/23.7/Linux_x86_64/23.7/compilers/bin/nvc++)
set(CMAKE_CXX_FLAGS "-mp=gpu -gpu=cc90 -Minfo=accel,opt,mp -Wall -DTWO_PROD_FMA")
set(CMAKE_BUILD_TYPE "Release")

project(QD-FFT CXX)

include_directories(qd)
include_directories(fft)
include_directories(bench_util)

add_subdirectory(bench)