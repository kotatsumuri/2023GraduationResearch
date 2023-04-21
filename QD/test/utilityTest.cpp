#include <cmath>
#include <gtest/gtest.h>
#include "utility.hpp"

TEST(quickTwoSum, zero) {
    double s, e;
    utility::quick_two_sum(0.0, 0.0, &s, &e);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e);
}

TEST(quickTwoSum, normal) {
    double s, e;
    double a = 1.0, b = std::pow(2, -52);
    std::cout << utility::to_reg_str(a) << std::endl;
    std::cout << utility::to_reg_str(b) << std::endl;
    utility::quick_two_sum(a, b, &s, &e);
    std::cout << utility::to_reg_str(s) << std::endl;   
    std::cout << utility::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);
    
    b = std::pow(2, -53);
    std::cout << utility::to_reg_str(a) << std::endl;
    std::cout << utility::to_reg_str(b) << std::endl;
    utility::quick_two_sum(a, b, &s, &e);
    std::cout << utility::to_reg_str(s) << std::endl;   
    std::cout << utility::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);
}