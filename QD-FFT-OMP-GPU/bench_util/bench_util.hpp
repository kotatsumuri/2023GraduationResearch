#pragma once
#include <omp.h>

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
    omp_set_num_threads(48);
#pragma omp parallel for reduction(+ : sum)
    for (uint64_t i = 0; i < n; i++) {
        sum += error_bit(real[i], actual[i]);
    }
    return sum / n;
}

double average_error_bit(uint64_t n, const qd_complex *real, const qd_complex *actual) {
    double sum = 0;
    omp_set_num_threads(48);
#pragma omp parallel for reduction(+ : sum)
    for (uint64_t i = 0; i < n; i++) {
        sum += error_bit(real[i].re, actual[i].re);
        sum += error_bit(real[i].im, actual[i].im);
    }
    return sum / (2 * n);
}

qd *warming_up(qd *a, qd *b, uint64_t n) {
    qd *c = (qd *)calloc(n, sizeof(qd));
#pragma omp target data map(to : a[0 : n], b[0 : n]) map(from : c[0 : n])
    {
#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            add(a[i], b[i], c[i]);
        }
    }
    return c;
}