#include <cmath>
#include <gtest/gtest.h>
#include <qd/dd_real.h>
#include "util.hpp"

TEST(quick_two_sum, zero) {
    double s, e;
    s = util::quick_two_sum(0.0, 0.0, &e);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e);
}

TEST(quick_two_sum, normal) {
    double s, e;
    double a = 1.0, b = std::pow(2, -52);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    s = util::quick_two_sum(a, b, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);
    
    b = std::pow(2, -53);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    s = util::quick_two_sum(a, b, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);
}

TEST(two_sum, zero) {
    double s, e;
    s = util::two_sum(0.0, 0.0, &e);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e);
}

TEST(two_sum, normal) {
    double s, e;
    double a = 1.0, b = std::pow(2, -52);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    s = util::two_sum(a, b, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);

    s = util::two_sum(b, a, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0 + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e);
    
    b = std::pow(2, -53);
    std::cout << util::to_reg_str(a) << std::endl;
    std::cout << util::to_reg_str(b) << std::endl;
    s = util::two_sum(a, b, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);

    s = util::two_sum(b, a, &e);
    std::cout << util::to_reg_str(s) << std::endl;   
    std::cout << util::to_reg_str(e) << std::endl;   
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e);
}

TEST(three_sum, zero) {
    double s, e0, e1;
    s = util::three_sum(0.0, 0.0, 0.0, &e0, &e1);
    EXPECT_EQ(0.0, s);
    EXPECT_EQ(0.0, e0);
    EXPECT_EQ(0.0, e1);
    s = util::three_sum(0.0, 0.0, 0.0, &e0);
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
    s = util::three_sum(a, b, c, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0 + std::pow(2, -52) + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e0);
    EXPECT_EQ(0.0, e1);
    s = util::three_sum(a, b, c, &e0);
    EXPECT_EQ(1.0 + std::pow(2, -52) + std::pow(2, -52), s);
    EXPECT_EQ(0.0, e0);
    
    c = std::pow(2, -53);
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    s = util::three_sum(a, b, c, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0 + std::pow(2, -51), s);
    EXPECT_EQ(- std::pow(2, -53), e0);
    EXPECT_EQ(0.0, e1);
    s = util::three_sum(a, b, c, &e0);
    EXPECT_EQ(1.0 + std::pow(2, -51), s);
    EXPECT_EQ(- std::pow(2, -53), e0);

    b = std::pow(2, -53);
    c = std::pow(2, -106);
    std::cout << "b:" << util::to_reg_str(b) << std::endl;
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    s = util::three_sum(a, b, c, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e0);
    EXPECT_EQ(std::pow(2, -106), e1);
    s = util::three_sum(a, b, c, &e0);
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53), e0);

    c = std::pow(2, -105);
    std::cout << "b:" << util::to_reg_str(b) << std::endl;
    std::cout << "c:" << util::to_reg_str(c) << std::endl;
    s = util::three_sum(a, b, c, &e0, &e1);
    std::cout << "s:" << util::to_reg_str(s) << std::endl;
    std::cout << "e0:" << util::to_reg_str(e0) << std::endl;
    std::cout << "e1:" << util::to_reg_str(e1) << std::endl;
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53) + std::pow(2, -105), e0);
    EXPECT_EQ(0.0, e1);
    s = util::three_sum(a, b, c, &e0);
    EXPECT_EQ(1.0, s);
    EXPECT_EQ(std::pow(2, -53) + std::pow(2, -105), e0);
}

TEST(dd_add_dd, normal) {
    dd_real a = dd_real::_pi;
    dd_real b = dd_real::_pi * dd_real(std::pow(2, -53));
    double s0, s1;
    util::dd_add_dd_dd(a.x[0], a.x[1], b.x[0], b.x[1], &s0, &s1);
    a += b;
    EXPECT_EQ(a.x[0], s0);
    EXPECT_EQ(a.x[1], s1);
    std::cout << util::to_reg_str(a.x[0]) << std::endl;
    std::cout << util::to_reg_str(s0) << std::endl;
    std::cout << util::to_reg_str(a.x[1]) << std::endl;
    std::cout << util::to_reg_str(s1) << std::endl;
}