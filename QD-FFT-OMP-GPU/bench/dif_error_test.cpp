#include <Arg.hpp>
#include <bench_util.hpp>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg       = Arg(argc, argv, 2, {"--debug", "--cos", "--range"}, {}, {});
    bool is_debug = arg.has_flag("--debug");
    bool is_cos   = arg.has_flag("--cos");
    bool is_range = arg.has_flag("--range");

    uint64_t p       = atoi(arg._argv[1]);
    uint64_t start_n = 1ull;
    uint64_t end_n   = 1ull << p;
    if (!is_range)
        start_n = end_n;

    std::cout << "n, real_ave_errorbit, imag_ave_errorbit" << std::endl;

    for (uint64_t n = start_n; n <= end_n; n <<= 1) {
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
        std::cout << n << ", ";
        std::cout << average_error_bit(n, actual_x, x) << ", ";
        std::cout << average_error_bit(n, actual_ix, ix) << std::endl;
    }
    return 0;
}