#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include "arg.hpp"
#include "bench_util.hpp"

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, 2, argv);
    bool is_debug = arg.has("--debug");
    bool is_cos   = arg.has("--cos");
    bool is_range = arg.has("--range");

    int p = atoi(arg.argv[1]);
    int n = 1ull << p;

    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table, cos_table);

    qd x[n];
    qd ix[n];
    if (is_cos) {
        copy_vector(n, cos_table, x);
        zero(n, ix);
    } else {
        rand_vector(n, x);
        rand_vector(n, ix);
    }
    qd actual_x[n];
    qd actual_ix[n];
    copy_vector(n, x, actual_x);
    copy_vector(n, ix, actual_ix);

    dif(n, x, ix, cos_table, sin_table);
    if (is_debug) {
        std::cout << "dif" << std::endl;
        std::cout << "x" << std::endl;
        print_vector(n, x);
        std::cout << "ix" << std::endl;
        print_vector(n, ix);
    }
    inv_dif(n, x, ix, cos_table, sin_table);
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