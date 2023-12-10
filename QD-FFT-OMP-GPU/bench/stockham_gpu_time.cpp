#include <omp.h>

#include <Arg.hpp>
#include <Timer.hpp>
#include <bench_util.hpp>
#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, argv, 2, {"--debug", "--range"}, {"-loops", "-start"}, {10, 1});
    bool is_debug = arg.has_flag("--debug");
    bool is_range = arg.has_flag("--range");

    uint64_t start_p = arg.get_option("-start");
    uint64_t end_p   = atoi(arg._argv[1]);
    uint64_t start_n = 1ull << start_p;
    uint64_t end_n   = 1ull << end_p;
    uint64_t K       = arg.get_option("-loops");

    if (!is_range) {
        start_p = end_p;
        start_n = end_n;
    }

    qd *base_w  = (qd *)calloc(end_n, sizeof(qd));
    qd *base_iw = (qd *)calloc(end_n, sizeof(qd));
    make_cos_table_gpu(end_n, base_w);
    make_sin_table(end_n, base_iw, base_w);

    if (warming_up(base_w, base_iw, 1 << 3) == NULL) {
        exit(1);
    } else {
        if (is_debug) {
            std::cout << "Warming up is done." << std::endl;
        }
    }

    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        Timer timer, h2d_timer, d2h_timer, kernel_timer;
        qd *x  = (qd *)calloc(n, sizeof(qd));
        qd *ix = (qd *)calloc(n, sizeof(qd));
        rand_vector(n, x);
        rand_vector(n, ix);

        qd *w  = base_w;
        qd *iw = base_iw;
        if (n != end_n) {
            w  = (qd *)calloc(n, sizeof(qd));
            iw = (qd *)calloc(n, sizeof(qd));
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++) {
                copy(base_w[end_n / n * i], w[i]);
                copy(base_iw[end_n / n * i], iw[i]);
            }
        }

        for (uint64_t k = 0; k < K; k++) {
            stockham_gpu(n, p, &x, &ix, w, iw, timer, h2d_timer, d2h_timer, kernel_timer);
        }

        std::cout << n << "," << timer.calc_ave_microsec() << ",";
        std::cout << h2d_timer.calc_ave_microsec() << ",";
        std::cout << d2h_timer.calc_ave_microsec() << ",";
        std::cout << kernel_timer.elapsed_microsec_once() / K << std::endl;
        free(x);
        free(ix);
        free(w);
        free(iw);
    }
    return 0;
}