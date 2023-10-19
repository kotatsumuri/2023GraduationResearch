#include "bench_util.hpp"
#include "fft.hpp"
#include "qd.hpp"

int main(int argc, char *argv[]) {
    if (argc < 2)
        return 1;
    int p = atoi(argv[1]);
    std::cout << "2^" << p << std::endl;
    int n = 1 << p;
    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table);
    qd x[n];
    qd ix[n];
    qd y[n];
    qd iy[n];
    qd correct[n];
    qd icorrect[n];
    int K = 100;

    double ave = 0;
    for (int t = 0; t < K; t++) {
        long long int error_bit_sum = 0;
        for (int i = 0; i < n; i++) {
            rand(x[i]);
            rand(ix[i]);
            copy(x[i], correct[i]);
            copy(ix[i], icorrect[i]);
        }

        stockham(n, p, x, ix, y, iy, cos_table, sin_table);

        if (!(p % 2))
            inv_stockham(n, p, x, ix, y, iy, cos_table, sin_table);
        else
            inv_stockham(n, p, y, iy, x, ix, cos_table, sin_table);

        for (int i = 0; i < n; i++) {
            error_bit_sum += error_bit(correct[i], x[i]);
            error_bit_sum += error_bit(icorrect[i], ix[i]);
        }
        ave += double(error_bit_sum) / double(2 * n);
    }

    std::cout << ave / K << std::endl;

    return 0;
}