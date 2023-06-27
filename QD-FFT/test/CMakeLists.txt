enable_testing()
include(GoogleTest)

add_executable(FFTTest fft_test.cpp)
target_link_libraries(FFTTest fft qd util_calc util_io)
target_compile_options(FFTTest PUBLIC -O2 -Wall)