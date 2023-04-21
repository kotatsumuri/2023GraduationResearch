enable_testing()

add_executable(UtilityTest utilityTest.cpp)
target_link_libraries(UtilityTest PRIVATE utility gtest_main)

include(GoogleTest)
gtest_discover_tests(UtilityTest)