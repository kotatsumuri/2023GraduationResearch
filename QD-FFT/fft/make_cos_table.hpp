#pragma once
#include <cstdint>
#include <queue>

#include "../qd/qd.hpp"

void make_quater_cos_table(uint32_t n, qd cos_table[]) {
    init(cos_table[0], 1.0);
    if (n <= 2)
        return;
    init(cos_table[n / 4], 0.0);
    if (n <= 4)
        return;
    sqrt(0.5, cos_table[n / 8]);
    if (n <= 8)
        return;
    std::queue<int> que;
    qd h;
    uint32_t i, j, k;
    que.push(n / 8);
    while (!que.empty()) {
        i = que.front();
        j = i / 2;
        k = n / 4 - j;
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

void make_cos_table(uint32_t n, qd cos_table[]) {
    qd qcos_table[n / 4 + 1];
    make_quater_cos_table(n, qcos_table);
    for(uint32_t i = 0;i < n / 4;i++) {
        copy(qcos_table[i], cos_table[i]);
        copy(qcos_table[n / 4 - i], cos_table[i + n / 4]);
        minus(cos_table[i + n / 4]);
        copy(qcos_table[i], cos_table[i + n / 2]);
        minus(cos_table[i + n / 2]);
        copy(qcos_table[n / 4 - i], cos_table[i + 3 * n / 4]);
    }
}

void make_sin_table(uint32_t n, qd sin_table[]) {
    qd qcos_table[n / 4 + 1];
    make_quater_cos_table(n, qcos_table);
    for(uint32_t i = 0;i < n / 4;i++) {
        copy(qcos_table[n / 4 - i], sin_table[i]);
        copy(qcos_table[i], sin_table[i + n / 4]);
        copy(qcos_table[n / 4 - i], sin_table[i + n / 2]);
        minus(sin_table[i + n / 2]);
        copy(qcos_table[i], sin_table[i + 3 * n / 4]);
        minus(sin_table[i + 3 * n / 4]);
    }
}