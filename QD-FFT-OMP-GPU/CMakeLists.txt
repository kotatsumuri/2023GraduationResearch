cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_COMPILER /system/apps/ubuntu/20.04-202304/nvhpc/23.7/Linux_x86_64/23.7/compilers/bin/nvc++)

set(CMAKE_CXX_FLAGS "-mp=gpu -gpu=cc80 -Minfo=accel,opt,mp -Wall -DTWO_PROD_FMA")

project(QD-FFT CXX)

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qd)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fft)

add_subdirectory(bench)