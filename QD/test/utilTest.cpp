#include <cmath>
#include <gtest/gtest.h>
#include "util.hpp"

TEST(quick_two_sum, zero) {
    double s, e;
    util::quick_two_sum(0.0, 0.0, &s, &e);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e);
}

TEST(quick_two_sum, normal) {
    double s, e;
    double a = 1.0, b = std::pow(2, -52);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    util::quick_two_sum(a, b, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);
    
    b = std::pow(2, -53);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    util::quick_two_sum(a, b, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);
}

TEST(two_sum, zero) {
    double s, e;
    util::two_sum(0.0, 0.0, &s, &e);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e);
}

TEST(two_sum, normal) {
    double s, e;
    double a = 1.0, b = std::pow(2, -52);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    util::two_sum(a, b, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);

    util::two_sum(b, a, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);
    
    b = std::pow(2, -53);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    util::two_sum(a, b, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);

    util::two_sum(b, a, &s, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);
}

TEST(three_sum, zero) {
    double s, e0, e1;
    util::three_sum(0.0, 0.0, 0.0, &s, &e0, &e1);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e0);
    EXPECT_EQ(0.0, e1);
    util::three_sum(0.0, 0.0, 0.0, &s, &e0);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e0);
}

TEST(three_sum, normal) {
    double a, b, c, s, e0, e1;
    a = 1.0;
    b = std::pow(2, -52);
    c = std::pow(2, -52);
    std::cout << "a:" << util::to_reg_str(a) << std::endl;
    std::cout << "b:" << util::to_reg_str(b) << std::endl;
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    util::three_sum(a, b, c, &s, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0 + std::pow(2, -52) + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e0);
    EXPECT_EQ(0.0, e1);
    util::three_sum(a, b, c, &s, &e0);
    EXPECT_EQ(1.0 + std::pow(2, -52) + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e0);
    
    c = std::pow(2, -53);
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    util::three_sum(a, b, c, &s, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0 + std::pow(2, -51), s);
    EXPECT_EQ(- std::pow(2, -53), e0);
    EXPECT_EQ(0.0, e1);
    util::three_sum(a, b, c, &s, &e0);
    EXPECT_EQ(1.0 + std::pow(2, -51), s);
    EXPECT_EQ(- std::pow(2, -53), e0);

    b = std::pow(2, -53);
    c = std::pow(2, -106);
    std::cout << "b:" << util::to_reg_str(b) << std::endl;
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    util::three_sum(a, b, c, &s, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e0);
    EXPECT_EQ(std::pow(2, -106), e1);
    util::three_sum(a, b, c, &s, &e0);
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e0);

    c = std::pow(2, -105);
    std::cout << "b:" << util::to_reg_str(b) << std::endl;
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    util::three_sum(a, b, c, &s, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53) + std::pow(2, -105), e0);
    EXPECT_EQ(0.0, e1);
    util::three_sum(a, b, c, &s, &e0);
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53) + std::pow(2, -105), e0);
}