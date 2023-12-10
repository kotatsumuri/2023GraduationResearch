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

    std::cout << "n, real-average-error-bit, imag-average-error-bit, average-error-bit" << std::endl;
    for (uint64_t n = start_n, p = start_p; n <= end_n; n <<= 1, p++) {
        qd *w = base_w;
        if (n != end_n) {
            w = (qd *)calloc(n / 4 + 1, sizeof(qd));
#pragma omp parallel for
            for (uint64_t i = 0; i < n / 4; i++) {
                copy(base_w[end_n / n * i], w[i]);
            }
        }
        qd *x  = (qd *)calloc(n, sizeof(qd));
        qd *ix = (qd *)calloc(n, sizeof(qd));
        rand_vector(n, x);
        rand_vector(n, ix);

        qd *actual_x  = (qd *)calloc(n, sizeof(qd));
        qd *actual_ix = (qd *)calloc(n, sizeof(qd));
        copy_vector(n, x, actual_x);
        copy_vector(n, ix, actual_ix);

        if (is_debug) {
            std::cout << "x" << std::endl;
            print_vector(n, x);
            std::cout << "ix" << std::endl;
            print_vector(n, ix);
        }

        stockham_gpu_n4(n, p, &x, &ix, w);
        if (is_debug) {
            std::cout << "stockham" << std::endl;
            std::cout << "x" << std::endl;
            print_vector(n, x);
            std::cout << "ix" << std::endl;
            print_vector(n, ix);
        }

        inv_stockham_gpu_n4(n, p, &x, &ix, w);
        if (is_debug) {
            std::cout << "inv_stockham" << std::endl;
            std::cout << "x" << std::endl;
            print_vector(n, x);
            std::cout << "ix" << std::endl;
            print_vector(n, ix);
        }

        double ave_error_bit_real = average_error_bit(n, actual_x, x);
        double ave_error_bit_imag = average_error_bit(n, actual_ix, ix);
        std::cout << n << ", ";
        std::cout << ave_error_bit_real << ", ";
        std::cout << ave_error_bit_imag << ", ";
        std::cout << (ave_error_bit_real + ave_error_bit_imag) / 2 << std::endl;

        free(x);
        free(ix);
        free(actual_x);
        free(actual_ix);
        free(w);
    }
    return 0;
}