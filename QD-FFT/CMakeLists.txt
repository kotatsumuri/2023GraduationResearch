cmake_minimum_required(VERSION 3.25)
project(QD-FFT CXX)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -DTWO_PROD_FMA")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qd)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fft)

add_subdirectory(bench)