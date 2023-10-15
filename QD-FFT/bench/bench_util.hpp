#pragma once
#include "qd.hpp"

double error_bit(qd real, qd actual) {
    qd c;
    sub(actual, real, c);
    fabs(c, c);
    int d            = (int)ceil(log2(c[0]));
    int log_ulp_real = log_ulp(real);
    return std::max(d - log_ulp_real + 1, 0);
}