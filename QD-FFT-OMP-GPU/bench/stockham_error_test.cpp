#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include "Arg.hpp"
#include "bench_util.hpp"

int main(int argc, char *argv[]) {
    Arg arg          = Arg(argc, 2, argv);
    bool is_debug    = arg.has("--debug");
    bool is_cos_real = arg.has("--cos-real");
    bool is_cos_imag = arg.has("--cos-imag");

    uint64_t p = atoi(arg.argv[1]);
    uint64_t n = 1ull << p;

    qd w[n], iw[n];
    make_cos_table(n, w);
    make_sin_table(n, iw, w);

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

    qd actual_x[n];
    qd actual_ix[n];
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
        std::cout << "dif" << std::endl;
        std::cout << "x" << std::endl;
        print_vector(n, x);
        std::cout << "ix" << std::endl;
        print_vector(n, ix);
    }

    inv_stockham(n, p, &x, &ix, w, iw);
    if (is_debug) {
        std::cout << "inv_dif" << std::endl;
        std::cout << "x" << std::endl;
        print_vector(n, x);
        std::cout << "ix" << std::endl;
        print_vector(n, ix);
    }

    std::cout << average_error_bit(n, actual_x, x) << std::endl;
    std::cout << average_error_bit(n, actual_ix, ix) << std::endl;

    return 0;
}