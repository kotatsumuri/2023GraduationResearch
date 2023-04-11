add_executable(testUtility.out testUtility.c)
target_link_libraries(testUtility.out Utility)

add_executable(testDFT.out testDFT.c)
target_link_libraries(testDFT.out Utility)
target_link_libraries(testDFT.out DFT)