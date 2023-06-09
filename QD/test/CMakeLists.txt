enable_testing()
include(GoogleTest)

add_executable(UtilTest utilTest.cpp)
target_compile_options(UtilTest PUBLIC -O2 -Wall)
target_link_libraries(UtilTest PRIVATE gtest_main qd util)
gtest_discover_tests(UtilTest)


add_executable(QDTest qdTest.cpp)
target_compile_options(QDTest PRIVATE -O2 -Wall)
target_link_libraries(QDTest PRIVATE gtest_main QD qd util)
gtest_discover_tests(QDTest)

add_executable(PrecisionTest precisionTest.cpp)
target_compile_options(PrecisionTest PRIVATE -O2 -Wall)
target_link_libraries(PrecisionTest PRIVATE QD qd util)
