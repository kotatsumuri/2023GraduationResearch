#include <iostream>

#include "fft.hpp"
#include "qd.hpp"

int main() {
    int p = 3;
    int n = 1 << p;
    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table);
    qd x[n];
    qd ix[n];
    qd y[n];
    qd iy[n];
    for(int i = 0;i < n;i++) {
        copy(cos_table[i], x[i]);
        copy(cos_table[i], ix[i]);
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(ix[i]) << std::endl;
    }
    // for(int i = 0;i < n;i++) {
    //     init(x[i], 1);
    //     init(ix[i], 1);
    // }
    stockham(n, p, x, ix, y, iy, cos_table, sin_table);
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(x[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(ix[i]) << std::endl;
    // }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(y[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(iy[i]) << std::endl;
    }

    if (!(p % 2))
        inv_stockham(n, p, x, ix, y, iy, cos_table, sin_table);
    else
        inv_stockham(n, p, y, iy, x, ix, cos_table, sin_table);

    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(cos_table[i]) << std::endl;
        std::cout << to_bin_string(x[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(cos_table[i]) << std::endl;
        std::cout << to_bin_string(ix[i]) << std::endl;
    }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(y[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(iy[i]) << std::endl;
    // }

    return 0;
}