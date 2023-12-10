#pragma once
#include "../qd/qd.hpp"

void make_quater_cos_table_gpu(uint64_t n, qd cos_table[]) {
    uint64_t n4 = n >> 2;
    uint64_t p  = (uint64_t)log2(n);
#pragma omp target data map(from : cos_table[0 : n4 + 1])
    {
#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i <= n / 4; i++) {
            double sign[32];
            if (i == 0) {
                init(cos_table[0], 1);
            } else if (i == n4) {
                zero(cos_table[i]);
            } else {
                uint64_t K = 0;
                uint64_t k = i;
                while (k != n4) {
                    k <<= 1;
                    sign[K] = 2;
                    if (k > n4) {
                        k       = (n >> 1) - k;
                        sign[K] = -2;
                    }
                    K++;
                }
                zero(cos_table[i]);
                qd tmp;
                for (uint64_t j = K; j > 0; j--) {
                    mul_pwr2(cos_table[i], sign[j - 1], tmp);
                    add(tmp, 2, tmp);
                    sqrt(tmp, cos_table[i]);
                    div_pwr2(cos_table[i], 2, cos_table[i]);
                }
            }
        }
    }
}

void make_cos_table_gpu(uint64_t n, qd cos_table[]) {
    if (n == 2) {
        init(cos_table[0], 1);
        init(cos_table[1], 0);
        return;
    }
    make_quater_cos_table_gpu(n, cos_table);
#pragma omp parallel for
    for (uint64_t i = 0; i < (n >> 2); i++) {
        copy(cos_table[(n >> 2) - i], cos_table[i + (n >> 2)]);
        minus(cos_table[i + (n >> 2)]);
        copy(cos_table[i], cos_table[i + (n >> 1)]);
        minus(cos_table[i + (n >> 1)]);
        copy(cos_table[(n >> 2) - i], cos_table[i + 3ull * (n >> 2)]);
    }
}