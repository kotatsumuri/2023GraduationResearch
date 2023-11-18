#pragma once
#include <qd.hpp>

double error_bit(const qd real, const qd actual) {
    qd c;
    sub(actual, real, c);
    fabs(c, c);
    int d            = (int)ceil(log2(c[0]));
    int log_ulp_real = log_ulp(real);
    return std::max(d - log_ulp_real + 1, 0);
}

double average_error_bit(uint64_t n, const qd *real, const qd *actual) {
    double sum = 0;
    for (uint64_t i = 0; i < n; i++) {
        sum += error_bit(real[i], actual[i]);
    }
    return sum / n;
}