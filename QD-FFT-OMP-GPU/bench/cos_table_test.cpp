#include <fft.hpp>
#include <iostream>
#include <qd.hpp>

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