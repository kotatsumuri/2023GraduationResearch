cmake_minimum_required(VERSION 3.25)

project(QD-FFT CXX)

set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -pg -g")

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_subdirectory(src)
add_subdirectory(test)