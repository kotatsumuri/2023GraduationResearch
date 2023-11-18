#pragma once
#include "qd.hpp"

void copy_vector(uint64_t n, const qd *x, qd *y) {
    for (uint64_t i = 0; i < n; i++) {
        copy(x[i], y[i]);
    }
}

void init_vector(uint64_t n, qd *x, double a) {
    for (uint64_t i = 0; i < n; i++) {
        init(x[i], a);
    }
}

void print_vector(uint64_t n, const qd *x) {
    for (uint64_t i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i]) << std::endl;
    }
}