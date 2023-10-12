#include <iostream>

#include "fft.hpp"
#include "qd.hpp"

int main() {
    int n = 8;
    int p = 3;
    qd cos_table[n + 1];
    make_quater_cos_table(n, cos_table);
    qd x[n];
    qd ix[n];
    qd y[n];
    qd iy[n];
    for (int i = 0; i < n; i++) {
        init(x[i], 1);
        zero(ix[i]);
    }
    stockham(n, p, x, ix, y, iy, cos_table);
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(ix[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(y[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(iy[i]) << std::endl;
    }

    if (!(p % 2))
        inv_stockham(n, p, x, ix, y, iy, cos_table);
    else
        inv_stockham(n, p, y, iy, x, ix, cos_table);

    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(x[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(ix[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(y[i]) << std::endl;
    }
    for (int i = 0; i < n; i++) {
        std::cout << to_bin_string(iy[i]) << std::endl;
    }

    return 0;
}