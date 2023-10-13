#include "bench_util.hpp"
#include "fft.hpp"
#include "qd.hpp"

int main() {
    int p = 3;
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

    double ave = 0;
    for (int t = 0; t < 1; t++) {
        long long int error_bit_sum = 0;

        for(int i = 0;i < n;i++) {
            copy(cos_table[i], x[i]);
            copy(cos_table[i], ix[i]);
        }
        // rand(x[0]);
        // rand(ix[0]);
        // // copy(x[0], ix[0]);
        // copy(x[0], correct[0]);
        // copy(ix[0], icorrect[0]);
        for (int i = 0; i < n; i++) {
            // rand(x[i]);
            // rand(ix[i]);
            // init(x[i], 1);
            // init(ix[i], 1);
            // copy(x[0], ix[i]);
            // copy(ix[0], x[i]);
            copy(x[i], correct[i]);
            copy(ix[i], icorrect[i]);
        }

        stockham(n, p, x, ix, y, iy, cos_table, sin_table);

        if (!(p % 2))
            inv_stockham(n, p, x, ix, y, iy, cos_table, sin_table);
        else
            inv_stockham(n, p, y, iy, x, ix, cos_table, sin_table);
        
        for(int i = 0;i < n;i++) {
            std::cout << to_bin_string(x[i]) << std::endl;
            std::cout << to_bin_string(correct[i]) << std::endl;
            std::cout << error_bit(correct[i], x[i]) << std::endl;
            error_bit_sum += error_bit(correct[i], x[i]);
            error_bit_sum += error_bit(icorrect[i], ix[i]);
        }
        ave += double(error_bit_sum) / double(2 * n);
    }

    std::cout << ave / 1 << std::endl;
    
    return 0;
}