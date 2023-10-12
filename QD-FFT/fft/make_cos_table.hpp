#pragma once
#include <cstdint>
#include <queue>

#include "../qd/qd.hpp"

void make_quater_cos_table(uint16_t n, qd cos_table[]) {
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
    uint16_t i, j, k;
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