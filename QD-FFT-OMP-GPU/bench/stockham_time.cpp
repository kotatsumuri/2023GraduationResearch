#include <omp.h>

#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include "Arg.hpp"
#include "Timer.hpp"

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, 2, argv);
    bool is_range = arg.has("--range");

    uint64_t start_p = 1ull;
    uint64_t end_p   = atoi(arg.argv[1]);
    uint64_t start_n = 1ull;
    uint64_t end_n   = 1ull << end_p;
    uint64_t K       = 10;

    if (!is_range) {
        start_p = end_p;
        start_n = end_n;
    }

    qd base_w[end_n];
    qd base_iw[end_n];
    make_cos_table(end_n, base_w);
    make_sin_table(end_n, base_iw, base_w);

    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        if (n != end_n) {
            qd w[n];
            qd iw[n];

            for (uint64_t i = 0; i < n; i++) {
                copy(base_w[end_n / n * i], w[i]);
                copy(base_iw[end_n / n * i], iw[i]);
            }
        }

        Timer timer;
        qd *x  = (qd *)calloc(n, sizeof(qd));
        qd *ix = (qd *)calloc(n, sizeof(qd));
        rand_vector(n, x);
        rand_vector(n, ix);

        if (n != end_n) {
            for (uint64_t k = 0; k < K; k++) {
                timer.start();
                stockham(n, p, &x, &ix, w, iw);
                timer.stop();
            }
        } else {
            for (uint64_t k = 0; k < K; k++) {
                timer.start();
                stockham(n, p, &x, &ix, base_w, base_iw);
                timer.stop();
            }
        }
        std::cout << n << "," << timer.calc_ave_microsec() << std::endl;
    }

    return 0;
}