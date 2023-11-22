#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

#include <bench_util.hpp>

int main(int argc, char *argv[]) {
    if (argc < 2)
        return 1;
    int p = atoi(argv[1]);
    std::cout << "2^" << p << std::endl;
    int n = 1 << p;
    double cos_table_n[n * 4];
    double sin_table_n[n * 4];
    uint64_t n1 = 1 << (p / 2);
    uint64_t n2 = 1 << ((p + 1) / 2);
    double cos_table_n1[n1 * 4];
    double sin_table_n1[n1 * 4];
    double cos_table_n2[n2 * 4];
    double sin_table_n2[n2 * 4];
    make_cos_table(n, cos_table_n);
    make_sin_table(n, sin_table_n, cos_table_n);
    for (uint64_t i = 0; i < n1; i++) {
        copy(cos_table_n + i * (n / n1) * 4, cos_table_n1 + i * 4);
        copy(sin_table_n + i * (n / n1) * 4, sin_table_n1 + i * 4);
    }
    for (uint64_t i = 0; i < n2; i++) {
        copy(cos_table_n + i * (n / n2) * 4, cos_table_n2 + i * 4);
        copy(sin_table_n + i * (n / n2) * 4, sin_table_n2 + i * 4);
    }
    double x[n * 8];
    double y[n * 8];
    double correct[n * 8];
    int K = 1;

    double ave = 0;
    for (int t = 0; t < K; t++) {
        long long int error_bit_sum = 0;
        for (int i = 0; i < n; i++) {
            rand(x + i * 8);
            rand(x + i * 8 + 4);
            copy(x + i * 8, correct + i * 8);
            copy(x + i * 8 + 4, correct + i * 8 + 4);
        }

        sixstep_fft(n, p, x, y, cos_table_n, sin_table_n, cos_table_n1, sin_table_n1, cos_table_n2, sin_table_n2);
        inv_stockham(n, p, y, x, cos_table_n, sin_table_n);

        for (int i = 0; i < n; i++) {
            error_bit_sum += error_bit(correct + i * 8, y + i * 8);
            error_bit_sum += error_bit(correct + i * 8 + 4, y + i * 8 + 4);
        }
        ave += double(error_bit_sum) / double(2 * n);
    }

    std::cout << ave / K << std::endl;

    return 0;
}