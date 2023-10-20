#include <omp.h>

#include <chrono>
#include <iostream>

#include "Timer.hpp"
#include "fft.hpp"
#include "qd.hpp"

int main(int argc, char *argv[]) {
    if (argc < 2)
        return 1;
    
    for(int p = 1;p <= atoi(argv[1]);p++) {
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
        for (int i = 0; i < 100; i++) {
            timer.start();
            stockham(n, p, x, ix, y, iy, cos_table, sin_table);
            timer.stop();
        }
        std::cout << n << "," << timer.calc_ave_microsec() << std::endl;
    }
    return 0;
}