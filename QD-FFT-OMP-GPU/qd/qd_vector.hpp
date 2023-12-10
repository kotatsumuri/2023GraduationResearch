#pragma once
#include <omp.h>

#include "qd.hpp"

void copy_vector(uint64_t n, const qd *x, qd *y) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        copy(x[i], y[i]);
    }
}

void copy_vector(uint64_t n, double *x[], qd *y) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        copy(x[i], y[i]);
    }
}

void copy_vector(uint64_t n, const qd_complex *x, qd_complex *y) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        copy(x[i].re, y[i].re);
        copy(x[i].im, y[i].im);
    }
}

void init_vector(uint64_t n, qd *x, double a) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        init(x[i], a);
    }
}

void print_vector(uint64_t n, const qd *x) {
    for (uint64_t i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i]) << std::endl;
    }
}

void print_vector(uint64_t n, const qd_complex *x) {
    for (uint64_t i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i].re) << std::endl;
        std::cout << to_bin_string(x[i].im) << std::endl;
    }
}

void rand_vector(uint64_t n, qd *x) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        rand(x[i]);
    }
}

void rand_vector(uint64_t n, qd_complex *x) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        rand(x[i].re);
        rand(x[i].im);
    }
}

void zero(uint64_t n, qd *x) {
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        zero(x[i]);
    }
}