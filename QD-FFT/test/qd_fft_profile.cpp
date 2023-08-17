#include "fft.hpp"
#include "qd.hpp"

int main() {
    const uint16_t N = (1 << 15);

    QD::qd cos_table[N / 4 + 1];
    QD::qd x[N];
    QD::qd ix[N];
    QD::qd y[N];
    QD::qd iy[N];

    FFT::make_cos_table(N, cos_table);

    for (uint16_t j = 0; j < N; j++) {
        QD::rand(x[j]);
        QD::rand(ix[j]);
    }

    FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);
}