add_executable(FFTTest fft_test.cpp)
target_link_libraries(FFTTest fft qd util_calc util_io)
target_compile_options(FFTTest PUBLIC -O2 -Wall)

add_executable(QDFFTPrecisionTest qd_fft_precision_test.cpp)
target_link_libraries(QDFFTPrecisionTest fft qd util_calc util_io)
target_compile_options(QDFFTPrecisionTest PUBLIC -O2 -Wall)

add_executable(QDFFTTime qd_fft_time.cpp)
target_link_libraries(QDFFTTime fft qd util_calc util_io)
target_compile_options(QDFFTTime PUBLIC -O2 -Wall)

add_executable(DoubleFFTTime double_fft_time.cpp)
target_link_libraries(DoubleFFTTime fft qd util_calc util_io)
target_compile_options(DoubleFFTTime PUBLIC -O2 -Wall)

add_executable(QDFFTErrorBits qd_fft_error_bits.cpp)
target_link_libraries(QDFFTErrorBits fft qd util_calc util_io)
target_compile_options(QDFFTErrorBits PUBLIC -O2 -Wall)

add_executable(QDFFTProfile qd_fft_profile.cpp)
target_link_libraries(QDFFTProfile fft qd util_calc util_io profiler tcmalloc)
target_compile_options(QDFFTProfile PUBLIC -O2 -Wall)