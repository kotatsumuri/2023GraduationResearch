#include <omp.h>

#include <chrono>
#include <iostream>

#include "Timer.hpp"
#include "fft.hpp"
#include "qd.hpp"

int main(int argc, char *argv[]) {
    if (argc < 3)
        return 1;
    int p = atoi(argv[1]);
    int n = 1 << p;
    Timer timer;
    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table);
    qd x[n];
    qd ix[n];
    qd y[n];
    qd iy[n];
    for (int i = 0; i < n; i++) {
        rand(x[i]);
        rand(ix[i]);
    }
    omp_set_num_threads(atoi(argv[2]));
    for (int i = 0; i < 100; i++) {
        timer.start();
        stockham(n, p, x, ix, y, iy, cos_table, sin_table);
        timer.stop();
    }
    std::cout << timer.calc_ave_msec() << std::endl;

    return 0;
}