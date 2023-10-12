#include <cmath>
#include <iostream>
#include <vector>

#include "fft.hpp"
#include "qd.hpp"
#include "util_calc.hpp"

int main() {
    const uint16_t M = 15;
    const uint16_t K = 10;
    uint16_t N       = 2;
    std::vector<int> bin(1025, 0);

    for (uint16_t i = 0; i < M; i++, N *= 2) {
        QD::qd cos_table[N / 4 + 1];
        FFT::make_cos_table(N, cos_table);

        for (uint16_t k = 0; k < K; k++) {
            QD::qd x[N];
            QD::qd ix[N];
            QD::qd x_[N];
            QD::qd ix_[N];
            QD::qd y[N];
            QD::qd iy[N];
            
            for (uint16_t j = 0; j < N; j++) {
                QD::rand(x[j]);
                QD::rand(ix[j]);
                QD::copy(x[j], x_[j]);
                QD::copy(ix[j], ix_[j]);
            }

            FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);
            FFT::inv_decimation_in_frequency(N, N, x, ix, y, iy, cos_table);

            for (long long int j = 0; j < N; j++) {
                QD::sub(x[j], x_[j], x_[j]);
                QD::sub(ix[j], ix_[j], ix_[j]);
                bin[-util::calc::ufp(x_[j][0])]++;
                bin[-util::calc::ufp(ix_[j][0])]++;
            }
        }
    }

    for(int i = 0;i < 1025;i++) {
        std::cout << -i << "," << bin[i] << std::endl;
    }
}