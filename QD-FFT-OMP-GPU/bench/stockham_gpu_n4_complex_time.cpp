#include <omp.h>

#include <Arg.hpp>
#include <Timer.hpp>
#include <bench_util.hpp>
#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, argv, 2, {"--debug", "--range"}, {"-loops", "-start"}, {10, 2});
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

    qd *base_w = (qd *)calloc(end_n / 4 + 1, sizeof(qd));
    make_quater_cos_table_gpu(end_n, base_w);

    if (warming_up(base_w, base_w, 1 << 3) == NULL) {
        exit(1);
    } else {
        if (is_debug) {
            std::cout << "Warming up is done." << std::endl;
        }
    }
    std::cout << "n, average-time, h2d-time, d2h-time, kernel-time, hd2-time(%), d2h-time(%), kernel-time(%)" << std::endl;
    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        Timer timer, h2d_timer, d2h_timer, kernel_timer;
        qd_complex *x = (qd_complex *)calloc(n, sizeof(qd_complex));
        rand_vector(n, x);

        qd *w = base_w;
        if (n != end_n) {
            w = (qd *)calloc(n / 4 + 1, sizeof(qd));
#pragma omp parallel for
            for (uint64_t i = 0; i < n / 4; i++) {
                copy(base_w[end_n / n * i], w[i]);
            }
        }

        for (uint64_t k = 0; k < K; k++) {
            stockham_gpu_n4_complex(n, p, &x, w, timer, h2d_timer, d2h_timer, kernel_timer);
        }

        std::cout << n << "," << timer.calc_ave_microsec() << ",";
        std::cout << h2d_timer.calc_ave_microsec() << ",";
        std::cout << d2h_timer.calc_ave_microsec() << ",";
        std::cout << kernel_timer.elapsed_microsec_once() / K << ",";
        std::cout << h2d_timer.calc_ave_microsec() / timer.calc_ave_microsec() * 100 << ",";
        std::cout << d2h_timer.calc_ave_microsec() / timer.calc_ave_microsec() * 100 << ",";
        std::cout << kernel_timer.elapsed_microsec_once() / K / timer.calc_ave_microsec() * 100 << std::endl;
        free(x);
        free(w);
    }
    return 0;
}