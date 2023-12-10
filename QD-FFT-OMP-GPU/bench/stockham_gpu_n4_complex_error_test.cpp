#include <Arg.hpp>
#include <bench_util.hpp>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, argv, 2, {"--debug", "--range"}, {}, {});
    bool is_debug = arg.has_flag("--debug");
    bool is_range = arg.has_flag("--range");

    uint64_t start_p = 2ull;
    uint64_t end_p   = atoi(arg._argv[1]);
    uint64_t start_n = 1ull << start_p;
    uint64_t end_n   = 1ull << end_p;
    if (!is_range) {
        start_p = end_p;
        start_n = end_n;
    }

    qd *base_w = (qd *)calloc(end_n / 4 + 1, sizeof(qd));
    make_quater_cos_table_gpu(end_n, base_w);

    if (warming_up(base_w, base_w, end_n / 4 + 1) == NULL) {
        exit(1);
    } else {
        if (is_debug) {
            std::cout << "Warming up is done." << std::endl;
        }
    }

    std::cout << "n, average-error-bit" << std::endl;
    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        qd *w = base_w;
        if (n != end_n) {
            w = (qd *)calloc(n / 4 + 1, sizeof(qd));
#pragma omp parallel for
            for (uint64_t i = 0; i < n / 4; i++) {
                copy(base_w[end_n / n * i], w[i]);
            }
        }
        qd_complex *x = (qd_complex *)calloc(n, sizeof(qd_complex));
        rand_vector(n, x);

        qd_complex *actual_x = (qd_complex *)calloc(n, sizeof(qd_complex));
        copy_vector(n, x, actual_x);

        if (is_debug) {
            std::cout << "x" << std::endl;
            print_vector(n, x);
        }

        stockham_gpu_n4_complex(n, p, &x, w);
        if (is_debug) {
            std::cout << "stockham" << std::endl;
            std::cout << "x" << std::endl;
            print_vector(n, x);
        }

        inv_stockham_gpu_n4_complex(n, p, &x, w);
        if (is_debug) {
            std::cout << "inv_stockham" << std::endl;
            std::cout << "x" << std::endl;
            print_vector(n, x);
        }

        double ave_error_bit_real = average_error_bit(n, actual_x, x);
        std::cout << n << ", ";
        std::cout << ave_error_bit_real << std::endl;

        free(x);
        free(actual_x);
        free(w);
    }
    return 0;
}