#include <Arg.hpp>
#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

int main(int argc, char *argv[]) {
    Arg arg(argc, argv, 2, {"--debug", "--range"}, {}, {});
    bool debug_flag = arg.has_flag("--debug");
    bool range_flag = arg.has_flag("--range");

    uint64_t p       = atoi(arg._argv[1]);
    uint64_t start_n = 1ull;
    uint64_t end_n   = 1ull << p;
    if (!range_flag)
        start_n = end_n;

    std::cout << "n, max_error, min_error, ave_error" << std::endl;
    for (uint64_t n = start_n; n <= end_n; n <<= 1) {
        qd *cos_table = (qd *)calloc(n / 4 + 1, sizeof(qd));
        make_quater_cos_table_gpu(n, cos_table);

        double max_error = 0;
        double min_error = 1e100;
        double ave_error = 0;
        for (uint64_t i = 0; i <= std::max(n >> 4, uint64_t(1)); i++) {
            qd cos2, sin2, one, error;
            sqr(cos_table[i], cos2);
            sqr(cos_table[n / 4 - i], sin2);
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
            print_vector(n / 4, cos_table);
            std::cout << std::endl;
        }
        std::cout << n << ", " << max_error << ", " << min_error << ", " << ave_error / (n >> 2) << std::endl;
        free(cos_table);
    }
    return 0;
}