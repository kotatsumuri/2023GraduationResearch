#pragma once
#include <iostream>

namespace QD_Lib {
    class QD {        
        public:
            double x[4];
            QD();
            QD(double x0);
            QD(double x0, double x1, double x2, double x3);
            void renormalize();
            void renormalize(double a);
            static QD qd_add_d_qd(const QD& a, double b);
            static QD qd_add_qd_qd(const QD& a, const QD& b);

            static QD qd_sub_d_qd(const QD& a, double b);
            static QD qd_sub_qd_qd(const QD& a, const QD& b);
            static QD d_sub_qd_qd(double a, const QD& b);

            static QD qd_mul_d_qd(const QD& a, double b);
            static QD qd_mul_qd_qd(const QD& a, const QD& b);
            static QD qd_mul_pwr2_qd(const QD& a, double b);
            static QD d_mul_d_qd(double a, double b);

            static QD qd_div_d_qd(const QD& a, double b);
            static QD qd_div_qd_qd(const QD& a, const QD& b);
            static QD d_div_qd_qd(double a, const QD& b);
            static QD d_div_d_qd(double a, double b);

            static QD minus(const QD& a);
            static QD pow(const QD& a, int n);
            static QD sqrt(const QD& a);
            static QD root(const QD& a, int n);
            static QD exp(const QD& a);
            static QD log(const QD& a);
            static QD cos(unsigned long long int k, unsigned long long int n);
            static QD sin(unsigned long long int k, unsigned long long int n);
            

            QD& operator =(const QD& r);
            QD& operator =(double l);
            QD operator +(const QD& r);
            QD operator +(double r);
            QD operator -(const QD& r);
            QD operator -(double r);
            const QD operator -();
            QD operator *(const QD& r);
            QD operator *(double r);
            QD operator /(const QD& r);
            QD operator /(double r);
            QD operator ^(int r);
            bool operator ==(const QD& r);
            bool operator !=(const QD& r);
            bool operator > (const QD& r);
            bool operator >=(const QD& r);
            bool operator < (const QD& r);
            bool operator <=(const QD& r);
            operator double();
    };
    QD operator+(double l, const QD& r);
    QD operator-(double l, const QD& r);
    QD operator*(double l, const QD& r);
    QD operator/(double l, const QD& r);
    std::ostream& operator <<(std::ostream& os, const QD& a);
}