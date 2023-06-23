#include "qd.hpp"

using namespace QD_Lib;

QD& QD::operator =(const QD& r) {
    x[0] = r.x[0];
    x[1] = r.x[1];
    x[2] = r.x[2];
    x[3] = r.x[3];
    return *this;
}

QD& QD::operator =(double r) {
    x[0] = r;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    return *this;
}

QD QD::operator +(const QD& r) {
    return qd_add_qd_qd(*this, r);
}

QD QD::operator +(double r) {
    return qd_add_d_qd(*this, r);
}

QD QD_Lib::operator +(double l, const QD& r) {
    return QD::qd_add_d_qd(r, l);
}

QD QD::operator -(const QD& r) {
    return qd_sub_qd_qd(*this, r);
}

QD QD::operator -(double r) {
    return qd_add_d_qd(*this, -r);
}

QD QD_Lib::operator -(double l, const QD& r) {
    return QD::d_sub_qd_qd(l, r);
}

const QD QD::operator -() {
    return QD(-x[0], -x[1], -x[2], -x[3]);
}

QD QD::operator *(const QD& r) {
    return qd_mul_qd_qd(*this, r);
}

QD QD::operator *(double r) {
    return qd_mul_d_qd(*this, r);
}

QD QD_Lib::operator *(double l, const QD& r) {
    return QD::qd_mul_d_qd(r, l);
}

QD QD::operator /(const QD& r) {
    return qd_div_qd_qd(*this, r);
}

QD QD::operator /(double r) {
    return qd_div_d_qd(*this, r);
}

QD QD_Lib::operator /(double l, const QD& r) {
    return QD::d_div_qd_qd(l, r);
}

QD QD::operator ^(int r) {
    return pow(*this, r);
}

bool QD::operator ==(const QD& r) {
    return (x[0] != r.x[0] ? 0 :
            x[1] != r.x[1] ? 0 :
            x[2] != r.x[2] ? 0 : 
            x[3] == r.x[3]);
}

bool QD::operator !=(const QD& r) {
    return (x[0] != r.x[0] ? 1 :
            x[1] != r.x[1] ? 1 :
            x[2] != r.x[2] ? 1 : 
            x[3] != r.x[3]);
}

bool QD::operator >(const QD& r) {
    return (x[0] != r.x[0] ? x[0] > r.x[0] :
            x[1] != r.x[1] ? x[1] > r.x[1] :
            x[2] != r.x[2] ? x[2] > r.x[2] : 
            x[3] > r.x[3]);
}

bool QD::operator >=(const QD& r) {
    return (x[0] != r.x[0] ? x[0] > r.x[0] :
            x[1] != r.x[1] ? x[1] > r.x[1] :
            x[2] != r.x[2] ? x[2] > r.x[2] : 
            x[3] >= r.x[3]);
}

bool QD::operator <(const QD& r) {
    return (x[0] != r.x[0] ? x[0] < r.x[0] :
            x[1] != r.x[1] ? x[1] < r.x[1] :
            x[2] != r.x[2] ? x[2] < r.x[2] : 
            x[3] < r.x[3]);
}

bool QD::operator <=(const QD& r) {
    return (x[0] != r.x[0] ? x[0] < r.x[0] :
            x[1] != r.x[1] ? x[1] < r.x[1] :
            x[2] != r.x[2] ? x[2] < r.x[2] : 
            x[3] <= r.x[3]);
}

QD::operator double() {
    return x[0];
}