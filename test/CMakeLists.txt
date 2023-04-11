add_executable(testUtility.out testUtility.c)
target_link_libraries(testUtility.out Utility)
target_compile_options(testUtility.out PUBLIC -O2 -Wall)

add_executable(testDFT.out testDFT.c)
target_link_libraries(testDFT.out Utility)
target_link_libraries(testDFT.out DFT)
target_compile_options(testDFT.out PUBLIC -O2 -Wall)

add_executable(testFFT.out testFFT.c)
target_link_libraries(testFFT.out Utility)
target_link_libraries(testFFT.out DFT)
target_link_libraries(testFFT.out FFT)
target_compile_options(testFFT.out PUBLIC -O2 -Wall)