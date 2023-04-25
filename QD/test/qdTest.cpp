#include <cmath>
#include <iostream>
#include <gtest/gtest.h>
#include <qd/qd_real.h>
#include "util.hpp"
#include "qd.hpp"

TEST(qd, renormalize) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));

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
}

TEST(qd, qd_prod_d_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    QD::qd_prod_d_qd(&a, std::pow(2, -53), &c);
    qd_real b = qd_real(1.0) + qd_real(std::pow(2, -53));
    b *= std::pow(2, -53);
    std::cout << c << std::endl;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(b.x[i], c.x[i]);
}

TEST(qd, qd_prod_qd_qd) {
    QD a(1.0, std::pow(2, -53), 0.0, 0.0);
    QD b(1.0, std::pow(2, -53), 0.0, 0.0);
    QD c;
    QD::qd_prod_qd_qd(&a, &b, &c);
    qd_real d = qd_real(1.0) + qd_real(std::pow(2, -53));
    qd_real e = qd_real(1.0) + qd_real(std::pow(2, -53));
    d = d * e;
    std::cout << c << std::endl;
    for(int i = 0;i < 4;i++)
        EXPECT_EQ(d.x[i], c.x[i]);
}