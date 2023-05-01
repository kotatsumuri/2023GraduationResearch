#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <bitset>
#include "util.hpp"

std::string util::to_reg_str(double a) {
    unsigned long long row_bit = *(unsigned long long *)&a;
    std::bitset<64> bin(row_bit);
    std::ostringstream ret;
    ret << ((bool)(row_bit >> 63) ? "-2^" : "+2^") << std::left << std::setw(5) << std::to_string(((int)(row_bit >> 52) & 2047) - 1023);
    return ret.str() + "*1." + bin.to_string().substr(12);
}

void util::split(double a, double* a_hi, double* a_lo) {
    *a_hi = (double)((1 << 27) + 1) * a;
    *a_hi = *a_hi - (*a_hi - a);
    *a_lo = a - *a_hi;
}

double util::quick_two_sum(double a, double b, double* e) {
    double s = a + b;
    *e = b - (s - a);
    return s;
}

double util::two_sum(double a, double b, double* e) {
    double s = a + b;
    double v = s - a;
    *e = (a - (s - v)) + (b - v);
    return s;
}

double util::three_sum(double a, double b, double c, double* e0, double* e1) {
    double s = two_sum(a, b, e0);
    s = two_sum(s, c, e1);
    *e0 = two_sum(*e0, *e1, e1);
    return s;
}

double util::three_sum(double a, double b, double c, double* e0) {
    double e1;
    double s = two_sum(a, b, e0);
    s = two_sum(s, c, &e1);
    *e0 += e1;
    return s;
}

double util::two_prod(double a, double b, double* e) {
    double a_hi, a_lo, b_hi, b_lo;
    double p = a * b;
    split(a, &a_hi, &a_lo);
    split(b, &b_hi, &b_lo);
    *e = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return p;
}

double util::six_three_sum(double a, double b, double c, double d, double e, double f, double* e0, double* e1) {
    double t[3];
    double s = three_sum(a, b, c, e0, e1);
    t[0] = three_sum(d, e, f, &t[1], &t[2]);
    s = two_sum(s, t[0], &t[0]);
    *e0 = two_sum(*e0, t[1], &t[1]);
    *e0 = two_sum(*e0, t[0], &t[0]);
    *e1 += t[2];
    *e1 += t[1];
    *e1 += t[0];
    return s;
}

void util::dd_add_dd_dd(double a0, double a1, double b0, double b1, double* s0, double* s1) {
    double t[2];
    *s0 = two_sum(a0, b0, s1);
    t[0] = two_sum(a1, b1, &t[1]);
    *s1 += t[0];
    *s0 = quick_two_sum(*s0, *s1, s1);
    *s1 += t[1];
    *s0 = quick_two_sum(*s0, *s1, s1);
}

double util::nine_two_sum(double a, double b, double c, double d, double e, double f, double g, double h, double i, double* e0) {
    double t[6];
    double s = util::two_sum(a, b, e0);
    t[0] = util::two_sum(c, d, &t[1]);
    dd_add_dd_dd(s, *e0, t[0], t[1], &s, e0);
    t[0] = two_sum(e, f, &t[1]);
    t[2] = two_sum(g, h, &t[3]);
    dd_add_dd_dd(t[0], t[1], t[2], t[3], &t[0], &t[1]);
    dd_add_dd_dd(s, *e0, t[0], t[1], &s, e0);
    return three_sum(i, s, *e0, e0);
}