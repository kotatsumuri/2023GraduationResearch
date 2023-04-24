cmake_minimum_required(VERSION 3.26)

project(QD CXX)

# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
)

set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

include_directories(/usr/local/include)
link_directories(/usr/local/lib) 

add_subdirectory(src)
add_subdirectory(test)