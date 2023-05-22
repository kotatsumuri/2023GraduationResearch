#pragma once
#include <iostream>

class QD {
    friend std::ostream& operator << (std::ostream& os, const QD& a);
    
    public:
        double x[4];
        QD();
        QD(double x0, double x1, double x2, double x3);
        void renormalize();
        void renormalize(double a);
        static QD qd_add_d_qd(const QD& a, double b);
        static QD qd_add_qd_qd(const QD& a, const QD& b);
        static QD qd_sub_qd_qd(const QD& a, const QD& b);
        static QD d_sub_qd_qd(double a, const QD& b);
        static QD qd_mul_d_qd(const QD& a, double b);
        static QD qd_mul_qd_qd(const QD& a, const QD& b);
        static QD d_mul_d_qd(double a, double b);
        static QD qd_div_qd_qd(const QD& a, const QD& b);
        static QD qd_div_d_qd(const QD& a, double b);
        static QD d_div_qd_qd(double a, const QD& b);
        static QD d_div_d_qd(double a, double b);
        static QD pow(const QD& a, int n);
        static QD sqrt(const QD& a);
        static QD root(const QD& a, int n);
        QD& operator =(const QD& r);
        QD operator +(double r);
        QD operator +(const QD& r);
        QD operator -(double r);
        QD operator -();
        QD operator *(double r);
        QD operator *(const QD& r);
        QD operator /(const QD& r);
        QD operator ^(int r);
};