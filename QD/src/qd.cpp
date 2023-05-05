#include <iostream>
#include "qd.hpp"
#include "util.hpp"

std::ostream& operator << (std::ostream& os, const QD& a) {
    os << util::to_reg_str(a.x[0]) << "|" << util::to_reg_str(a.x[1]) << "|" << util::to_reg_str(a.x[2]) << "|" << util::to_reg_str(a.x[3]);
    return os;
}

QD::QD() {
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
}
QD::QD(double x0, double x1, double x2, double x3) {
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    x[3] = x3;
    renormalize();
}

void QD::renormalize() {
    double s, e, t[5];
    int k, i;

    t[4] = 0.0;
    s = util::quick_two_sum(x[2], x[3], &t[3]);
    s = util::quick_two_sum(x[1], s, &t[2]);
    t[0] = util::quick_two_sum(x[0], s, &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        s = util::quick_two_sum(s, t[i], &e);
        if(e != 0) {
            x[k++] = s;
            s = e;
        }
    }

    if(k < 4) {
        x[k++] = s;
        while(k < 4)
            x[k++] = 0.0;
    }
}

void QD::renormalize(double a) {
    double s, e, t[5];
    int k, i;

    s = util::quick_two_sum(x[3], a, &t[4]);
    s = util::quick_two_sum(x[2], s, &t[3]);
    s = util::quick_two_sum(x[1], s, &t[2]);
    t[0] = util::quick_two_sum(x[0], s, &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        s = util::quick_two_sum(s, t[i], &e);
        if(e != 0) {
            x[k++] = s;
            s = e;
        }
    }

    if(k < 4) {
        x[k++] = s;
        while(k < 4)
            x[k++] = 0.0;
    }
}

QD QD::qd_add_d_qd(const QD& a, double b) {
    double e;
    QD c;
    c.x[0] = util::two_sum(a.x[0], b, &e);
    c.x[1] = util::two_sum(a.x[1], e, &e);
    c.x[2] = util::two_sum(a.x[2], e, &e);
    c.x[3] = util::two_sum(a.x[3], e, &e);

    c.renormalize(e);
    return c;
}

QD QD::qd_add_qd_qd(const QD& a, const QD& b) {
    double t[4];
    QD c;
    c.x[0] = util::two_sum(a.x[0], b.x[0], &t[0]);
    c.x[1] = util::two_sum(a.x[1], b.x[1], &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_sum(a.x[2],  b.x[2], &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = util::two_sum(a.x[3],  b.x[3], &t[3]);
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[1] += t[0] + t[3];
    c.renormalize(t[1]);
    return c;
}

QD QD::qd_sub_qd_qd(const QD& a, const QD& b) {
    double t[4];
    QD c;
    c.x[0] = util::two_sum(a.x[0], -b.x[0], &t[0]);
    c.x[1] = util::two_sum(a.x[1], -b.x[1], &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_sum(a.x[2],  -b.x[2], &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = util::two_sum(a.x[3],  -b.x[3], &t[3]);
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[1] += t[0] + t[3];
    c.renormalize(t[1]);
    return c;
}

QD QD::qd_mul_d_qd(const QD& a, double b) {
    double t[4];
    QD c;
    c.x[0] = util::two_prod(a.x[0], b, &t[0]);
    c.x[1] = util::two_prod(a.x[1], b, &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_prod(a.x[2], b, &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = a.x[3] * b;
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[0] += t[1];
    c.renormalize(t[0]);
    return c;
}

QD QD::qd_mul_qd_qd(const QD& a, const QD& b) {
    double t[17];
    QD c;
    c.x[0] = util::two_prod(a.x[0], b.x[0], &c.x[1]); // 0, 1  
    t[0] = util::two_prod(a.x[0], b.x[1], &t[1]); // 1, 2 
    t[2] = util::two_prod(a.x[1], b.x[0], &t[3]); // 1, 2 
    c.x[1] = util::three_sum(c.x[1], t[0], t[2], &t[0], &t[2]); // 1, 2, 3 
    c.x[2] = util::two_prod(a.x[0], b.x[2], &t[4]); // 2, 3 
    t[5] = util::two_prod(a.x[1], b.x[1], &t[6]); // 2, 3
    t[7] = util::two_prod(a.x[2], b.x[0], &t[8]); // 2, 3 
    c.x[2] = util::six_three_sum(t[1], t[3], t[0], c.x[2], t[5], t[7], &t[0], &t[1]); // 2, 3, 4 
    c.x[3] = util::two_prod(a.x[0], b.x[3], &t[3]); // 3, 4
    t[5] = util::two_prod(a.x[1], b.x[2], &t[7]);  // 3, 4
    t[9] = util::two_prod(a.x[2], b.x[1], &t[10]);  // 3, 4
    t[11] = util::two_prod(a.x[3], b.x[0], &t[12]);  // 3, 4
    c.x[3] = util::nine_two_sum(t[2], t[4], t[6], t[8], t[0], c.x[3], t[5], t[9], t[11], &t[0]); // 3, 4 
    t[0] = a.x[1] * b.x[3] + a.x[2] * b.x[2] + a.x[3] * b.x[1] + t[1] + t[3] + t[7] + t[10] + t[12] + t[0];
    c.renormalize(t[0]);
    return c;
}

QD QD::qd_div_qd_qd(const QD& a, const QD& b) {
    double q;
    QD t0, t1;
    QD c;
    c.x[0] = a.x[0] / b.x[0];
    t0 = qd_mul_d_qd(b, c.x[0]);
    t0 = qd_sub_qd_qd(a, t0);
    c.x[1] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[1]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[2] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[2]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[3] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[3]);
    t0 = qd_sub_qd_qd(t0, t1);
    q = t0.x[0] / b.x[0];
    c.renormalize(q);
    return c;
}

QD QD::pow(const QD& a, int n) {
    QD c = a;
    QD b = QD(1.0, 0, 0, 0);
    while(n != 0) {
        if(n & 1)
            b = qd_mul_qd_qd(b, c);
        c = qd_mul_qd_qd(c, c);
        n >>= 1;
    }
    return b;
}

QD& QD::operator =(const QD& r) {
    x[0] = r.x[0];
    x[1] = r.x[1];
    x[2] = r.x[2];
    x[3] = r.x[3];
    return *this;
}

QD QD::operator +(double r) {
    return qd_add_d_qd(*this, r);
}

QD QD::operator +(const QD& r) {
    return qd_add_qd_qd(*this, r);
} 

QD QD::operator -(double r) {
    return qd_add_d_qd(*this, -r);
}

QD QD::operator -() {
    return QD(-x[0], -x[1], -x[2], -x[3]);
}

QD QD::operator *(double r) {
    return qd_mul_d_qd(*this, r);
}

QD QD::operator *(const QD& r) {
    return qd_mul_qd_qd(*this, r);
}

QD QD::operator /(const QD& r) {
    return qd_div_qd_qd(*this, r);
}

QD QD::operator ^(int r) {
    return pow(*this, r);
}