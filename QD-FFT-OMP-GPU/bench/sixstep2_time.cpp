#include <omp.h>

#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include "Timer.hpp"

int main(int argc, char *argv[]) {
    if (argc < 2)
        return 1;

    double base_cos_table[(1 << atoi(argv[1])) * 4];
    double base_sin_table[(1 << atoi(argv[1])) * 4];
    make_cos_table(1 << atoi(argv[1]), base_cos_table);
    make_sin_table(1 << atoi(argv[1]), base_sin_table, base_cos_table);

    for (int p = 1; p <= atoi(argv[1]); p++) {
        int n = 1 << p;
        Timer timer;
        double cos_table_n[n * 4];
        double sin_table_n[n * 4];
        uint64_t n1 = 1 << (p / 2);
        uint64_t n2 = 1 << ((p + 1) / 2);
        double cos_table_n1[n1 * 4];
        double sin_table_n1[n1 * 4];
        double cos_table_n2[n2 * 4];
        double sin_table_n2[n2 * 4];
        for (uint64_t i = 0; i < n; i++) {
            copy(base_cos_table + i * (1 << atoi(argv[1])) / n * 4, cos_table_n + i * 4);
            copy(base_sin_table + i * (1 << atoi(argv[1])) / n * 4, sin_table_n + i * 4);
        }
        for (uint64_t i = 0; i < n1; i++) {
            copy(cos_table_n + i * (n / n1) * 4, cos_table_n1 + i * 4);
            copy(sin_table_n + i * (n / n1) * 4, sin_table_n1 + i * 4);
        }
        for (uint64_t i = 0; i < n2; i++) {
            copy(cos_table_n + i * (n / n2) * 4, cos_table_n2 + i * 4);
            copy(sin_table_n + i * (n / n2) * 4, sin_table_n2 + i * 4);
        }
        double x[n * 8];
        double y[n * 8];

        for (int i = 0; i < n; i++) {
            rand(x + i * 8);
            rand(x + i * 8 + 4);
        }
        for (int i = 0; i < 10; i++) {
            timer.start();
            sixstep_fft(n, p, x, y, cos_table_n, sin_table_n, cos_table_n1, sin_table_n1, cos_table_n2, sin_table_n2);
            timer.stop();
        }
        std::cout << n << "," << timer.calc_ave_microsec() << std::endl;
    }
    return 0;
}