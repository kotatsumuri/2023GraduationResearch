#include <cmath>
#include <iostream>
#include <gtest/gtest.h>
#include <qd/qd_real.h>
#include "util.hpp"
#include "qd.hpp"

TEST(qd, renormalize) {
    QD a(1.0, std::pow(2, -52), 0.0, 0.0);
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -52));
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
    a = QD(M_PI, M_PI, M_PI, M_PI);
    b = M_PI * 4;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
}

TEST(qd, qd_add_d_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    a = QD::qd_add_d_qd(a, std::pow(2, -53));
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));
    b += std::pow(2, -53);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
    b = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = b.x[i];
    a = QD::qd_add_d_qd(a, M_PI);
    b += M_PI;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
}

TEST(qd, qd_add_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(std::pow(2, -53), 0, 0, 0);
    QD c;
    c = QD::qd_add_qd_qd(a, b);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    d += qd_real(std::pow(2, -53));
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
    d = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = d.x[i];
    }
    c = QD::qd_add_qd_qd(a, b);
    d += qd_real::_pi;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
}

TEST(qd, qd_sub_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(std::pow(2, -53), 0, 0, 0);
    QD c;
    c = QD::qd_sub_qd_qd(a, b);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    d -= qd_real(std::pow(2, -53));
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
    d = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = d.x[i];
    }
    c = QD::qd_sub_qd_qd(a, b);
    d -= qd_real::_pi;
    for(int i = 0;i < 4;i++) 
        EXPECT_EQ(d.x[i], c.x[i]);
}


TEST(qd, qd_prod_d_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    c = QD::qd_mul_d_qd(a, std::pow(2, -53));
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));
    b *= std::pow(2, -53);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], c.x[i]);
    b = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = b.x[i];
    b *= M_PI;
    c = QD::qd_mul_d_qd(a, M_PI);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], c.x[i]);
}

TEST(qd, qd_prod_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    c = QD::qd_mul_qd_qd(a, b);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    qd_real e = qd_real(1.0) + qd_real(std::pow(2, -53));
    d = qd_real::accurate_mul(d, e);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);

    d = qd_real::_pi;
    e = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = e.x[i];
    }
    d = qd_real::accurate_mul(d, e);
    c = QD::qd_mul_qd_qd(a, b);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
}

TEST(qd, qd_div_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    c = QD::qd_div_qd_qd(a, b);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    qd_real e = qd_real(1.0) + qd_real(std::pow(2, -53));
    d = qd_real::accurate_div(d, e);
    for(int i = 0;i < 4;i++) 
        EXPECT_EQ(d.x[i], c.x[i]);

    d = qd_real::_pi;
    e = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = e.x[i];
    }
    d = qd_real::accurate_div(d, e);
    c = QD::qd_div_qd_qd(a, b);
    for(int i = 0;i < 4;i++) 
        EXPECT_EQ(d.x[i], c.x[i]);
}

TEST(qd, pow) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b;
    b = QD::pow(a, 10);
    qd_real c = qd_real(1.0) + qd_real(std::pow(2, -53));   
    c = c^10;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(c.x[i], b.x[i]);
    c = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = c.x[i];
    c = c^100;
    b = QD::pow(a, 100);   
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(c.x[i], b.x[i]);
}

TEST(qd, sqrt) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b = QD::sqrt(a);
    qd_real c = qd_real(1.0) + qd_real(std::pow(2, -53));
    c = QD_API::sqrt(c);
    for(int i = 0;i < 4;i++) {
        EXPECT_EQ(c.x[i], b.x[i]);
        std::cout << util::to_reg_str(c.x[i]) << std::endl;
        std::cout << util::to_reg_str(b.x[i]) << std::endl;
    }
    c = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = c.x[i];
    c = QD_API::sqrt(c);
    b = QD::sqrt(a);
    for(int i = 0;i < 4;i++) {
        EXPECT_EQ(c.x[i], b.x[i]);
        std::cout << util::to_reg_str(c.x[i]) << std::endl;
        std::cout << util::to_reg_str(b.x[i]) << std::endl;
    }
}

TEST(qd, root) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b = QD::root(a, 3);
    qd_real c = qd_real(1.0) + qd_real(std::pow(2, -53));
    c = QD_API::nroot(c, 3);
    for(int i = 0;i < 4;i++) {
        EXPECT_EQ(c.x[i], b.x[i]);
        std::cout << util::to_reg_str(c.x[i]) << std::endl;
        std::cout << util::to_reg_str(b.x[i]) << std::endl;
    }
    c = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = c.x[i];
    c = QD_API::nroot(c, 3);
    b = QD::root(a, 3);
    for(int i = 0;i < 4;i++) {
        EXPECT_EQ(c.x[i], b.x[i]);
        std::cout << util::to_reg_str(c.x[i]) << std::endl;
        std::cout << util::to_reg_str(b.x[i]) << std::endl;
    }

    c = 8;
    for(int i = 0;i < 4;i++)
        a.x[i] = c.x[i];
    c = QD_API::nroot(c, 3);
    b = QD::root(a, 3);
    for(int i = 0;i < 4;i++) {
        EXPECT_EQ(c.x[i], b.x[i]);
        std::cout << util::to_reg_str(c.x[i]) << std::endl;
        std::cout << util::to_reg_str(b.x[i]) << std::endl;
    }
}