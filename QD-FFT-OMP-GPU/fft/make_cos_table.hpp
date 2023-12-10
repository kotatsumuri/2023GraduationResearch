#pragma once
#include <omp.h>

#include <cstdint>
#include <queue>

#include "../qd/qd.hpp"

void make_quater_cos_table(uint64_t n, qd cos_table[]) {
    init(cos_table[0], 1.0);
    if (n <= 2)
        return;
    init(cos_table[n >> 2], 0.0);
    if (n <= 4)
        return;
    sqrt(0.5, cos_table[n >> 3]);
    if (n <= 8)
        return;
    std::queue<uint64_t> que;
    qd h;
    uint64_t i, j, k;
    que.push(n >> 3);
    while (!que.empty()) {
        i = que.front();
        j = i >> 1;
        k = (n >> 2) - j;
        que.pop();
        mul_pwr2(cos_table[i], 2, h);
        add(h, 2.0, cos_table[k]);
        sqrt(cos_table[k], cos_table[j]);
        mul_pwr2(cos_table[j], 0.5, cos_table[j]);
        sub(2.0, h, cos_table[k]);
        sqrt(cos_table[k], h);
        mul_pwr2(h, 0.5, cos_table[k]);
        if (!(j & 1)) {
            que.push(j);
            que.push(k);
        }
    }
}

void make_cos_table(uint64_t n, qd cos_table[]) {
    if (n == 2) {
        init(cos_table[0], 1);
        init(cos_table[1], 0);
        return;
    }
    make_quater_cos_table(n, cos_table);
    for (uint64_t i = 0; i < (n >> 2); i++) {
        copy(cos_table[(n >> 2) - i], cos_table[i + (n >> 2)]);
        minus(cos_table[i + (n >> 2)]);
        copy(cos_table[i], cos_table[i + (n >> 1)]);
        minus(cos_table[i + (n >> 1)]);
        copy(cos_table[(n >> 2) - i], cos_table[i + 3ull * (n >> 2)]);
    }
}

void make_sin_table(uint64_t n, qd sin_table[]) {
    qd qcos_table[n / 4 + 1];
    if (n == 2) {
        init(sin_table[0], 0);
        init(sin_table[1], 0);
        return;
    }
    make_quater_cos_table(n, qcos_table);
    for (uint64_t i = 0; i < n / 4; i++) {
        copy(qcos_table[n / 4 - i], sin_table[i]);
        copy(qcos_table[i], sin_table[i + n / 4]);
        copy(qcos_table[n / 4 - i], sin_table[i + n / 2]);
        minus(sin_table[i + n / 2]);
        copy(qcos_table[i], sin_table[i + 3 * n / 4]);
        minus(sin_table[i + 3 * n / 4]);
    }
}

void make_sin_table(uint64_t n, qd sin_table[], qd cos_table[]) {
    if (n == 2) {
        init(sin_table[0], 0);
        init(sin_table[1], 0);
        return;
    }
    omp_set_num_threads(48);
#pragma omp parallel for
    for (uint64_t i = 0; i < (n >> 2); i++) {
        copy(cos_table[(n >> 2) - i], sin_table[i]);
        copy(cos_table[i], sin_table[i + (n >> 2)]);
        copy(cos_table[(n >> 2) - i], sin_table[i + (n >> 1)]);
        minus(sin_table[i + (n >> 1)]);
        copy(cos_table[i], sin_table[i + 3ull * (n >> 2)]);
        minus(sin_table[i + 3ull * (n >> 2)]);
    }
}

void make_quater_cos_table(uint64_t n, double cos_table[]) {
    init(cos_table, 1.0);
    if (n <= 2)
        return;
    init(cos_table + (n / 4) * 4, 0.0);
    if (n <= 4)
        return;
    sqrt(0.5, cos_table + (n / 8) * 4);
    if (n <= 8)
        return;
    std::queue<int> que;
    qd h;
    uint64_t i, j, k;
    que.push(n / 8);
    while (!que.empty()) {
        i = que.front();
        j = i / 2;
        k = n / 4 - j;
        que.pop();
        mul_pwr2(cos_table + i * 4, 2, h);
        add(h, 2.0, cos_table + k * 4);
        sqrt(cos_table + k * 4, cos_table + j * 4);
        mul_pwr2(cos_table + j * 4, 0.5, cos_table + j * 4);
        sub(2.0, h, cos_table + k * 4);
        sqrt(cos_table + k * 4, h);
        mul_pwr2(h, 0.5, cos_table + k * 4);
        if (!(j & 1)) {
            que.push(j);
            que.push(k);
        }
    }
}

void make_cos_table(uint64_t n, double cos_table[]) {
    double qcos_table[(n / 4 + 1) * 4];
    if (n == 2) {
        init(cos_table, 1);
        init(cos_table + 4, 0);
        return;
    }
    make_quater_cos_table(n, qcos_table);
    for (uint64_t i = 0; i < n / 4; i++) {
        copy(qcos_table + i * 4, cos_table + i * 4);
        copy(qcos_table + (n / 4 - i) * 4, cos_table + (i + n / 4) * 4);
        minus(cos_table + (i + n / 4) * 4);
        copy(qcos_table + i * 4, cos_table + (i + n / 2) * 4);
        minus(cos_table + (i + n / 2) * 4);
        copy(qcos_table + (n / 4 - i) * 4, cos_table + (i + 3 * n / 4) * 4);
    }
}

void make_sin_table(uint64_t n, double sin_table[], double cos_table[]) {
    if (n == 2) {
        init(sin_table, 0);
        init(sin_table + 4, 0);
        return;
    }
    for (uint64_t i = 0; i < n / 4; i++) {
        copy(cos_table + (n / 4 - i) * 4, sin_table + i * 4);
        copy(cos_table + i * 4, sin_table + (i + n / 4) * 4);
        copy(cos_table + (n / 4 - i) * 4, sin_table + (i + n / 2) * 4);
        minus(sin_table + (i + n / 2) * 4);
        copy(cos_table + i * 4, sin_table + (i + 3 * n / 4) * 4);
        minus(sin_table + (i + 3 * n / 4) * 4);
    }
}