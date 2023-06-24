#include "util_calc.hpp"

#include <bitset>

namespace util::calc {
void split(double a, double* a_hi, double* a_lo) {
    *a_hi = (double)((1 << 27) + 1) * a;
    *a_hi = *a_hi - (*a_hi - a);
    *a_lo = a - *a_hi;
}

void quick_two_sum(double a, double b, double* s, double* e) {
    *s = a + b;
    *e = b - (*s - a);
}

void two_sum(double a, double b, double* s, double* e) {
    *s       = a + b;
    double v = *s - a;
    *e       = (a - (*s - v)) + (b - v);
}

void three_sum(double a, double b, double c, double* s, double* e0,
               double* e1) {
    two_sum(a, b, s, e0);
    two_sum(*s, c, s, e1);
    two_sum(*e0, *e1, e0, e1);
}

void three_sum(double a, double b, double c, double* s, double* e0) {
    double e1;
    two_sum(a, b, s, e0);
    two_sum(*s, c, s, &e1);
    *e0 += e1;
}

void two_prod(double a, double b, double* p, double* e) {
    double a_hi, a_lo, b_hi, b_lo;
    *p = a * b;
    split(a, &a_hi, &a_lo);
    split(b, &b_hi, &b_lo);
    *e = ((a_hi * b_hi - *p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
}

void six_three_sum(double a, double b, double c, double d, double e, double f,
                   double* s, double* e0, double* e1) {
    double t[3];
    three_sum(a, b, c, s, e0, e1);
    three_sum(d, e, f, t, t + 1, t + 2);
    two_sum(*s, t[0], s, t);
    two_sum(*e0, t[1], e0, t + 1);
    two_sum(*e0, t[0], e0, t);
    *e1 += t[2];
    *e1 += t[1];
    *e1 += t[0];
}

void dd_add_dd_dd(double a0, double a1, double b0, double b1, double* s0,
                  double* s1) {
    double t[2];
    two_sum(a0, b0, s0, s1);
    two_sum(a1, b1, t, t + 1);
    *s1 += t[0];
    quick_two_sum(*s0, *s1, s0, s1);
    *s1 += t[1];
    quick_two_sum(*s0, *s1, s0, s1);
}

void nine_two_sum(double a, double b, double c, double d, double e, double f,
                  double g, double h, double i, double* s, double* e0) {
    double t[6];
    two_sum(a, b, s, e0);
    two_sum(c, d, t, t + 1);
    dd_add_dd_dd(*s, *e0, t[0], t[1], s, e0);
    two_sum(e, f, t, t + 1);
    two_sum(g, h, t + 2, t + 3);
    dd_add_dd_dd(t[0], t[1], t[2], t[3], &t[0], &t[1]);
    dd_add_dd_dd(*s, *e0, t[0], t[1], s, e0);
    three_sum(i, *s, *e0, s, e0);
}

int ufp(double a) {
    unsigned long long row_bit = *(unsigned long long*)&a;
    std::bitset<64> bin(row_bit);
    return ((int)(row_bit >> 52) & 2047) - 1023;
}

int ulp(double a) { return ufp(a) - 52; }
}  // namespace util::calc