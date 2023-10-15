#include "fft.hpp"
#include "qd.hpp"

int main() {
    int p = 3;
    int n = 1 << p;
    qd cos_table[n];
    qd sin_table[n];
    make_cos_table(n, cos_table);
    make_sin_table(n, sin_table);

    qd sqrt_harf, harf;
    sqrt(0.5, sqrt_harf);
    sqr(sqrt_harf, harf);
    // std::cout << to_bin_string(harf) << std::endl;

    for (int i = 0; i < n / 4; i++) {
        qd cos2, sin2, one, error;
        sqr(cos_table[i], cos2);
        sqr(sin_table[i], sin2);
        add(cos2, sin2, one);
        sub(1, one, error);
        fabs(error, error);
        std::cout << i << " " << to_bin_string(error) << std::endl;
    }

    return 0;
}