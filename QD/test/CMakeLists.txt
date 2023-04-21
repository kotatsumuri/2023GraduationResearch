enable_testing()

add_executable(UtilTest utilTest.cpp)
target_link_libraries(UtilTest PRIVATE util gtest_main)

include(GoogleTest)
gtest_discover_tests(UtilTest)