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
    ret << ((bool)(row_bit >> 63) ? "-2^" : "+2^") << std::setw(5) << std::to_string(((int)(row_bit >> 52) & 2047) - 1023);
    return ret.str() + "*1." + bin.to_string().substr(12);
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

void util::three_sum(double a, double b, double c, double* s, double* e) {
    double e1;
    two_sum(a,  b, s, e);
    two_sum(*s, c, s, &e1);
    *e += e1;
}