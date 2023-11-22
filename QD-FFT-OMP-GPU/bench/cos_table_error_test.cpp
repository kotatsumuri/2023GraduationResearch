#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include <Arg.hpp>

int main(int argc, char *argv[]) {
    Arg arg(argc, 2, argv);

    bool debug_flag = arg.has("--debug");
    bool range_flag = arg.has("--range");

    uint64_t p       = atoi(arg.argv[1]);
    uint64_t start_n = 1ull;
    uint64_t end_n   = 1ull << p;
    if (!range_flag)
        start_n = end_n;

    for (uint64_t n = start_n; n <= end_n; n <<= 1) {
        std::cout << "n: " << n << std::endl;
        qd cos_table[n];
        qd sin_table[n];
        make_cos_table(n, cos_table);
        make_sin_table(n, sin_table, cos_table);

        double max_error = 0;
        double min_error = 1e100;
        double ave_error = 0;
        for (uint64_t i = 0; i < std::max(n >> 2, uint64_t(1)); i++) {
            qd cos2, sin2, one, error;
            sqr(cos_table[i], cos2);
            sqr(sin_table[i], sin2);
            add(cos2, sin2, one);
            sub(1.0, one, error);
            fabs(error, error);
            max_error = std::max(max_error, error[0]);
            min_error = std::min(min_error, error[0]);
            ave_error += error[0];

            if (debug_flag)
                std::cout << "error: " << to_bin_string(error) << std::endl;
        }

        if (debug_flag) {
            std::cout << "cos_table:" << std::endl;
            print_vector(n, cos_table);
            std::cout << std::endl;

            std::cout << "sin_table:" << std::endl;
            print_vector(n, sin_table);
            std::cout << std::endl;
        }

        std::cout << "max_error: " << max_error << std::endl;
        std::cout << "min_error: " << min_error << std::endl;
        std::cout << "ave_error: " << ave_error / (n >> 2) << std::endl;
    }
    return 0;
}