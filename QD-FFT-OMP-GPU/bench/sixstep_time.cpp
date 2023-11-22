#include <omp.h>

#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include <Timer.hpp>

int main(int argc, char *argv[]) {
    if (argc < 2)
        return 1;

    uint64_t max_n = 1ull << atoi(argv[1]);
    qd base_cos_table[max_n];
    qd base_sin_table[max_n];
    make_cos_table(max_n, base_cos_table);
    make_sin_table(max_n, base_sin_table, base_cos_table);

    for (uint64_t p = 1ull; p <= atoi(argv[1]); p++) {
        uint64_t n = 1ull << p;
        Timer timer;
        qd cos_table[n];
        qd sin_table[n];
        uint64_t u = 1ull << (atoi(argv[1]) - p);
        for (uint64_t i = 0; i < n; i++) {
            uint64_t idx = i * u;
            copy(base_cos_table[idx], cos_table[i]);
            copy(base_sin_table[idx], sin_table[i]);
        }

        qd x[n];
        qd ix[n];
        qd y[n];
        qd iy[n];
        for (uint64_t i = 0; i < n; i++) {
            rand(x[i]);
            rand(ix[i]);
        }

        for (int i = 0; i < 10; i++) {
            timer.start();
            sixstep_fft(n, p, x, ix, y, iy, cos_table, sin_table);
            timer.stop();
        }
        std::cout << n << "," << timer.calc_ave_microsec() << std::endl;
    }
    return 0;
}