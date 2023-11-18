#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

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
    for (int i = 0; i < n; i++) {
        copy(cos_table[i], x[i]);
        copy(cos_table[i], ix[i]);
    }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(x[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(ix[i]) << std::endl;
    // }

    // for (int i = 0; i < n; i++) {
    //     // init(x[i], 1);
    //     init(ix[i], 0);
    // }

    double *x_[n];
    double *ix_[n];
    double *y_[n];
    double *iy_[n];

    for (int i = 0; i < n; i++) {
        x_[i]  = x[i];
        ix_[i] = ix[i];
        y_[i]  = y[i];
        iy_[i] = iy[i];
    }
    dif(n, n, x_, ix_, y_, iy_, cos_table, sin_table);

    // for (int i = 0; i < n; i++) {
    //     // std::cout << "cos_table[" << i << "] " << to_bin_string(cos_table[i]) << std::endl;
    //     std::cout << "        y[" << i << "] " << to_bin_string(x_[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(x[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(ix[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(y[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(iy[i]) << std::endl;
    // }

    inv_dif(n, n, x_, ix_, y_, iy_, cos_table, sin_table);
    for (int i = 0; i < n; i++) {
        std::cout << "cos_table[" << i << "] " << to_bin_string(cos_table[i]) << std::endl;
        std::cout << "        y[" << i << "] " << to_bin_string(y_[i]) << std::endl;
    }
    // for (int i = 0; i < n; i++) {
    //     // std::cout << to_bin_string(cos_table[i]) << std::endl;
    //     std::cout << "       iy[" << i << "] " << to_bin_string(iy[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(y[i]) << std::endl;
    // }
    // for (int i = 0; i < n; i++) {
    //     std::cout << to_bin_string(iy[i]) << std::endl;
    // }

    return 0;
}