#include <omp.h>

#include <Arg.hpp>
#include <Timer.hpp>
#include <fft.hpp>
#include <iostream>

#define THRESHOLD 10
#define THREADS 48

int main(int argc, char *argv[]) {
    Arg arg = Arg(argc, argv, 2, {}, {"-start"}, {1});
    Timer global_timer;
    uint64_t start_p = arg.get_option("-start");
    uint64_t end_p   = atoi(arg._argv[1]);
    uint64_t start_n = 1ull << start_p;
    uint64_t end_n   = 1ull << end_p;

    qd *base_w  = (qd *)calloc(end_n, sizeof(qd));
    qd *base_iw = (qd *)calloc(end_n, sizeof(qd));
    global_timer.start();
    make_cos_table(end_n, base_w);
    make_sin_table(end_n, base_iw, base_w);

    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        qd *x  = (qd *)calloc(n, sizeof(qd));
        qd *ix = (qd *)calloc(n, sizeof(qd));
        zero(n, x);
        zero(n, ix);

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

        omp_set_num_threads(THREADS);
        stockham_omp(n, p, &x, &ix, w, iw);

        for (uint32_t i = 1; i < (1 << 30); i <<= 1) {
            double min_elapsed_time = 1e10;
            for (uint32_t j = 0; j < 8; j++) {
                Timer timer;
                omp_set_num_threads(THREADS);
                for (uint32_t k = 0; k < i; k++) {
                    stockham_omp(n, p, &x, &ix, w, iw, timer);
                }
                if (timer.elapsed_time() < min_elapsed_time) {
                    min_elapsed_time = timer.elapsed_time();
                }
                if (timer.elapsed_time() < THRESHOLD)
                    break;
            }
            if (min_elapsed_time > THRESHOLD) {
                std::cout << n << ", " << min_elapsed_time / i;
                std::cout << std::endl;
                break;
            }
        }
        free(x);
        free(ix);
        free(w);
        free(iw);
    }
}