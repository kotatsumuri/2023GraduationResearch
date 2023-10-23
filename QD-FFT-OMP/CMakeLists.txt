cmake_minimum_required(VERSION 3.25)
project(QD-FFT CXX)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -DTWO_PROD_FMA -O2")

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qd)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fft)

add_subdirectory(bench)