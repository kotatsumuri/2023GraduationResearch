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

    copy_vector(n, cos_table, x);
    copy_vector(n, cos_table, ix);

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
    inv_dif(n, n, x_, ix_, y_, iy_, cos_table, sin_table);

    print_vector(n, x);
    print_vector(n, ix);

    return 0;
}