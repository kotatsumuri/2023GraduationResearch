#include <cmath>
#include <iostream>
#include <vector>

#include "fft.hpp"
#include "qd.hpp"
#include "util_calc.hpp"

int main() {
    const uint16_t M = 15;
    const uint16_t K = 10;
    uint16_t N       = 8;
    std::vector<int> bin(1025, 0);

    QD::qd sqrt_harf, harf;
    QD::sqrt(0.5, sqrt_harf);

    QD::qd cos_table[N / 4 + 1];
    FFT::make_cos_table(N, cos_table);

    QD::qd x[N];
    QD::qd ix[N];
    QD::qd x_[N];
    QD::qd ix_[N];
    QD::qd y[N];
    QD::qd iy[N];

    // for (uint16_t j = 0; j < N; j++) {
    //     QD::rand(x[j]);
    //     QD::rand(ix[j]);
    //     QD::copy(x[j], x_[j]);
    //     QD::copy(ix[j], ix_[j]);
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

    for (uint16_t j = 0; j < N; j++) {
        // QD::rand(x[j]);
        // QD::rand(ix[j]);
        QD::copy(x[j], x_[j]);
        QD::copy(ix[j], ix_[j]);
    }

    FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);
    FFT::inv_decimation_in_frequency(N, N, x, ix, y, iy, cos_table);

    for(int i = 0;i < N;i++) {
        std::cout << QD::to_bin_string(x[i]) << std::endl;
        std::cout << QD::to_bin_string(x_[i]) << std::endl;
    }

    // for (int i = 0; i < N / 4; i++) {
    //     QD::qd cos2, sin2, one, error;
    //     QD::sqr(cos_table[i], cos2);
    //     QD::sqr(cos_table[N / 4 - i], sin2);
    //     QD::add(cos2, sin2, one);
    //     QD::sub(1, one, error);
    //     QD::abs(error, error);
    //     std::cout << i << " " << QD::to_bin_string(error) << std::endl;
    // }

}