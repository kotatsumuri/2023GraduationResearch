#include "fft.hpp"

#include <cmath>
#include <iostream>

#include "qd.hpp"

const int N = 8;

int main() {
    QD::qd cos_table[N / 4 + 1];
    FFT::make_cos_table(N, cos_table);
    std::cout << "cos_table" << std::endl;
    for (int i = 0; i <= N / 4; i++) {
        std::cout << QD::to_bin_string(cos_table[i]) << std::endl;
    }

    QD::qd x[N], ix[N];
    QD::qd y[N], iy[N];

    QD::init(x[0], 1.0);
    QD::init(ix[0], 0.0);
    // QD::init(x[1], -1.0);
    // QD::init(ix[1], 0.0);
    // for (int i = 0; i < N; i++) {
    //     QD::init(x[i], 1.0);
    //     QD::init(ix[i], 0.0);
    // }
    QD::init(ix[0], 0);
    QD::init(ix[N / 4], 0);
    QD::init(ix[N / 2], 0);
    QD::init(ix[3 * N / 4], 0);
    QD::init(x[0], 1);
    QD::init(x[N / 4], 0);
    QD::init(x[N / 2], -1);
    QD::init(x[3 * N / 4], 0);

    for (int i = 1; i < N / 4; i++) {
        QD::copy(cos_table[i], x[i]);
        QD::copy(cos_table[i], x[N - i]);
        QD::init(x[N / 2 - i], -cos_table[i][0], -cos_table[i][1],
                 -cos_table[i][2], -cos_table[i][3]);
        QD::copy(x[N / 2 - i], x[N / 2 + i]);
        QD::init(ix[i], 0);
        QD::init(ix[N - i], 0);
        QD::init(ix[N / 2 + i], 0);
        QD::init(ix[N / 2 - i], 0);
    }

    std::cout << "x" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << QD::to_bin_string(x[i]) << std::endl;
        std::cout << QD::to_bin_string(ix[i]) << std::endl;
    }

    FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);

    std::cout << "FFT(x)" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << QD::to_bin_string(x[i]) << std::endl;
        std::cout << QD::to_bin_string(ix[i]) << std::endl;
    }

    FFT::inv_decimation_in_frequency(N, N, x, ix, y, iy, cos_table);

    std::cout << "iFFT(FFT(x))" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << QD::to_bin_string(x[i]) << std::endl;
        std::cout << QD::to_bin_string(ix[i]) << std::endl;
    }

    return 0;
}