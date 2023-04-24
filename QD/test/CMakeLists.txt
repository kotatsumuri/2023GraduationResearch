enable_testing()
include(GoogleTest)

add_executable(UtilTest utilTest.cpp)
target_link_libraries(UtilTest PRIVATE util gtest_main)
gtest_discover_tests(UtilTest)

add_executable(QDTest qdTest.cpp)
target_link_libraries(QDTest PRIVATE util gtest_main)
target_link_libraries(QDTest PRIVATE QD gtest_main)
target_link_libraries(QDTest PRIVATE qd gtest_main)
gtest_discover_tests(QDTest)