#include <fft.hpp>
#include <iostream>
#include <qd.hpp>
#include "bench_util.hpp"

int main(int argc, char *argv[]) {
    if(argc < 2)
        return 1;

    int p = atoi(argv[1]);
    int n = 1 << p;

    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table);
    
    qd x[n];
    qd ix[n];
    qd y[n];
    qd iy[n];

    rand_vector(n, x);
    rand_vector(n, ix);
    rand_vector(n, y);
    rand_vector(n, iy);

    qd actual_x[n];
    qd actual_ix[n];
    qd actual_y[n];
    qd actual_iy[n];

    copy_vector(n, x, actual_x);
    copy_vector(n, ix, actual_ix);
    copy_vector(n, y, actual_y);
    copy_vector(n, iy, actual_iy);

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

    average_error_bit(n, actual_x, x);
    average_error_bit(n, actual_ix, ix);
    average_error_bit(n, actual_y, y);
    average_error_bit(n, actual_iy, iy);

    return 0;
}