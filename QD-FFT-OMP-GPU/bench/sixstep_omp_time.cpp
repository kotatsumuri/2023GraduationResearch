#include <omp.h>

#include <Arg.hpp>
#include <Timer.hpp>
#include <chrono>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg(argc, argv, 4, {"--debug", "--range"}, {"-loops", "-start"}, {10, 1});
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
    make_cos_table(end_n, base_w);
    make_sin_table(end_n, base_iw, base_w);

    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        qd *w  = base_w;
        qd *iw = base_iw;
        if (n != end_n) {
            w  = (qd *)calloc(n, sizeof(qd));
            iw = (qd *)calloc(n, sizeof(qd));

            for (uint64_t i = 0; i < n; i++) {
                copy(base_w[end_n / n * i], w[i]);
                copy(base_iw[end_n / n * i], iw[i]);
            }
        }
        qd *x  = (qd *)calloc(n, sizeof(qd));
        qd *ix = (qd *)calloc(n, sizeof(qd));
        rand_vector(n, x);
        rand_vector(n, ix);
        for (int thread_idx = 0; thread_idx < atoi(arg._argv[2]); thread_idx++) {
            int thread_num = atoi(arg._argv[3 + thread_idx]);
            Timer timer;
            omp_set_num_threads(thread_num);
            for (uint64_t k = 0; k < K; k++) {
                sixstep_omp(n, p, &x, &ix, w, iw, timer);
            }
            std::cout << thread_num << ", " << n << "," << timer.calc_ave_microsec() << std::endl;
        }
        free(x);
        free(ix);
        free(w);
        free(iw);
    }

    return 0;
}