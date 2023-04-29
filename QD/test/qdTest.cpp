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
    QD::qd_add_d_qd(&a, std::pow(2, -53), &a);
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));
    b += std::pow(2, -53);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
    b = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = b.x[i];
    QD::qd_add_d_qd(&a, M_PI, &a);
    b += M_PI;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], a.x[i]);
}

TEST(qd, qd_add_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(std::pow(2, -53), 0, 0, 0);
    QD c;
    QD::qd_add_qd_qd(&a, &b, &c);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    d += qd_real(std::pow(2, -53));
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
    d = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = d.x[i];
    }
    QD::qd_add_qd_qd(&a, &b, &c);
    d += qd_real::_pi;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
}

TEST(qd, qd_sub_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(std::pow(2, -53), 0, 0, 0);
    QD c;
    QD::qd_sub_qd_qd(&a, &b, &c);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    d -= qd_real(std::pow(2, -53));
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
    d = qd_real::_pi;
    for(int i = 0;i < 4;i++) {
        a.x[i] = d.x[i];
        b.x[i] = d.x[i];
    }
    QD::qd_sub_qd_qd(&a, &b, &c);
    d -= qd_real::_pi;
    for(int i = 0;i < 4;i++) 
        EXPECT_EQ(d.x[i], c.x[i]);
}


TEST(qd, qd_prod_d_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    QD::qd_mul_d_qd(&a, std::pow(2, -53), &c);
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));
    b *= std::pow(2, -53);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], c.x[i]);
    b = qd_real::_pi;
    for(int i = 0;i < 4;i++)
        a.x[i] = b.x[i];
    b *= M_PI;
    QD::qd_mul_d_qd(&a, M_PI, &c);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], c.x[i]);
}

TEST(qd, qd_prod_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    QD::qd_mul_qd_qd(&a, &b, &c);
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
    QD::qd_mul_qd_qd(&a, &b, &c);
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
}

TEST(qd, qd_div_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    QD::qd_div_qd_qd(&a, &b, &c);
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
    QD::qd_div_qd_qd(&a, &b, &c);
    for(int i = 0;i < 4;i++) 
        EXPECT_EQ(d.x[i], c.x[i]);
}