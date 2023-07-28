#include <cmath>
#include <iostream>

#include "fft.hpp"
#include "qd.hpp"

int main() {
    long long int N = 2;
    for (int i = 0; i < 5; i++) {
        QD::qd cos_table[N / 4 + 1];
        QD::qd x_[N];
        QD::qd ix_[N];
        QD::qd x[N];
        QD::qd ix[N];
        QD::qd y[N];
        QD::qd iy[N];

        FFT::make_cos_table(N, cos_table);

        // QD::init(x[0], 1.0);
        // QD::init(ix[0], 0.0);
        // QD::init(ix[0], 0);
        // QD::init(ix[N / 4], 0);
        // QD::init(ix[N / 2], 0);
        // QD::init(ix[3 * N / 4], 0);
        // QD::init(x[0], 1);
        // QD::init(x[N / 4], 0);
        // QD::init(x[N / 2], -1);
        // QD::init(x[3 * N / 4], 0);

        // for (long long int j = 1; j < N / 4; j++) {
        //     QD::copy(cos_table[j], x[j]);
        //     QD::copy(cos_table[j], x[N - j]);
        //     QD::init(x[N / 2 - j], -cos_table[j][0], -cos_table[j][1],
        //              -cos_table[j][2], -cos_table[j][3]);
        //     QD::copy(x[N / 2 - j], x[N / 2 + j]);
        //     QD::init(ix[j], 0);
        //     QD::init(ix[N - j], 0);
        //     QD::init(ix[N / 2 + j], 0);
        //     QD::init(ix[N / 2 - j], 0);
        // }

        for (long long int j = 0; j < N; j++) {
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
            std::cout << N << " " << j << " " << QD::to_bin_string(x_[j])
                      << std::endl;
            std::cout << N << " " << j << " " << QD::to_bin_string(ix_[j])
                      << std::endl;
        }
        N *= 2;
    }
}