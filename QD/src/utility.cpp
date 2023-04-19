#include "utility.hpp"

/**
 * @brief Calculate sum and error.
 * 
 * @param a 
 * @param b 
 * @param s sum
 * @param e error
 */
void utility::two_sum(double a, double b, double* s, double* e) {
    *s = a + b;
    double v = *s - a;
    *e = (a - (*s - v)) + (b - v);
}