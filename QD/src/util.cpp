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

/**
 * @brief Calculate sum and error. (|a| >= |b|)
 * 
 * @param a 
 * @param b (|a| >= |b|)
 * @param s sum
 * @param e error
 */
void util::quick_two_sum(double a, double b, double* s, double* e) {
    *s = a + b;
    *e = b - (*s - a);
}

/**
 * @brief Calculate sum and error.
 * 
 * @param a 
 * @param b 
 * @param s sum
 * @param e error
 */
void util::two_sum(double a, double b, double* s, double* e) {
    *s = a + b;
    double v = *s - a;
    *e = (a - (*s - v)) + (b - v);
}

void util::three_sum(double a, double b, double c, double* s, double* e0, double* e1) {
    two_sum(a,   b,   s,  e0);
    two_sum(*s,  c,   s,  e1);
    two_sum(*e0, *e1, e0, e1);
}

void util::three_sum(double a, double b, double c, double* s, double* e0) {
    double e1;
    two_sum(a,  b, s, e0);
    two_sum(*s, c, s, &e1);
    *e0 += e1;
}

void util::two_prod(double a, double b, double* p, double* e) {
    double a_hi, a_lo, b_hi, b_lo;
    *p = a * b;
    split(a, &a_hi, &a_lo);
    split(b, &b_hi, &b_lo);
    *e = ((a_hi * b_hi - *p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
}

void util::six_three_sum(double a, double b, double c, double d, double e, double f, double* s, double* e0, double* e1) {
    double t[3];
    three_sum(a, b, c, s, e0, e1);
    three_sum(d, e, f, &t[0], &t[1], &t[2]);
    two_sum(*s, t[0], s, &t[0]);
    two_sum(*e0, t[1], e0, &t[1]);
    two_sum(*e0, t[0], e0, &t[0]);
    *e1 += t[2];
    *e1 += t[1];
    *e1 += t[0];
}

void util::dd_add_dd_dd(double a0, double a1, double b0, double b1, double* s0, double* s1) {
    double t[2];
    two_sum(a0, b0, s0, s1);
    two_sum(a1, b1, &t[0], &t[1]);
    *s1 += t[0];
    quick_two_sum(*s0, *s1, s0, s1);
    *s1 += t[1];
    quick_two_sum(*s0, *s1, s0, s1);
}

void util::nine_two_sum(double a, double b, double c, double d, double e, double f, double g, double h, double i, double* s, double* e0) {
    double t[6];
    util::two_sum(a, b, s, e0);
    util::two_sum(c, d, &t[0], &t[1]);
    util::dd_add_dd_dd(*s, *e0, t[0], t[1], s, e0);
    util::two_sum(e, f, &t[0], &t[1]);
    util::two_sum(g, h, &t[2], &t[3]);
    util::dd_add_dd_dd(t[0], t[1], t[2], t[3], &t[0], &t[1]);
    util::dd_add_dd_dd(*s, *e0, t[0], t[1], s, e0);
    util::three_sum(i, *s, *e0, s, e0);
}