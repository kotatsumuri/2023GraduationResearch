#include <Arg.hpp>
#include <bench_util.hpp>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg          = Arg(argc, argv, 2, {"--debug", "--cos-real", "--cos-imag", "--range"}, {}, {});
    bool is_debug    = arg.has_flag("--debug");
    bool is_cos_real = arg.has_flag("--cos-real");
    bool is_cos_imag = arg.has_flag("--cos-imag");
    bool is_range    = arg.has_flag("--range");

    uint64_t start_p = 1ull;
    uint64_t end_p   = atoi(arg._argv[1]);
    uint64_t start_n = 1ull << start_p;
    uint64_t end_n   = 1ull << end_p;
    if (!is_range) {
        start_p = end_p;
        start_n = end_n;
    }

    qd *base_w  = (qd *)calloc(end_n, sizeof(qd));
    qd *base_iw = (qd *)calloc(end_n, sizeof(qd));
    make_cos_table(end_n, base_w);
    make_sin_table(end_n, base_iw, base_w);

    std::cout << "n, real-average-error-bit, imag-average-error-bit, average-error-bit" << std::endl;
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
        zero(n, x);
        zero(n, ix);
        if (is_cos_real) {
            copy_vector(n, w, x);
        }
        if (is_cos_imag) {
            copy_vector(n, w, ix);
        }
        if (!is_cos_real && !is_cos_imag) {
            rand_vector(n, x);
            rand_vector(n, ix);
        }

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

        stockham(n, p, &x, &ix, w, iw);
        if (is_debug) {
            std::cout << "stockham" << std::endl;
            std::cout << "x" << std::endl;
            print_vector(n, x);
            std::cout << "ix" << std::endl;
            print_vector(n, ix);
        }

        inv_stockham(n, p, &x, &ix, w, iw);
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
    }
    return 0;
}