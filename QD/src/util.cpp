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