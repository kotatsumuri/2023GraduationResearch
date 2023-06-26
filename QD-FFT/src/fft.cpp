#include "fft.hpp"

#include "qd.hpp"

namespace FFT {

void decimal_in_frequency(int initN, int n, QD::qd x[], QD::qd ix[], QD::qd y[],
                          QD::qd iy[], QD::qd cos_table[]) {
    if (n <= 1)
        return;

    const int harf_n = n / 2;
    int i            = 1;
    int j            = 1 + harf_n;
    // y[0] = x[0] + x[n / 2]
    QD::add(x[0], x[harf_n], y[0]);
    QD::add(ix[0], ix[harf_n], iy[0]);
    // y[n / 2] = x[0] - x[n / 2]
    QD::sub(x[0], x[harf_n], y[harf_n]);
    QD::sub(ix[0], ix[harf_n], iy[harf_n]);

    if (n > 2) {
        int i_                = harf_n - 1;
        int j_                = n - 1;
        const int harf_harf_n = harf_n / 2;
        const int initN_n     = initN / n;
        const int initN_4     = initN / 4;

        for (; i < harf_harf_n;) {
            QD::add(x[i], x[j], y[i]);
            QD::add(ix[i], ix[j], iy[i]);

            QD::sub(x[i], x[j], y[j]);
            QD::sub(ix[i], ix[j], iy[j]);

            double* a = (double*)cos_table[i * initN_n];
            double* b = (double*)cos_table[initN_4 - i * initN_n];

            QD::mul(y[j], a, x[i]);
            QD::mul(y[j], b, ix[i]);
            QD::mul(iy[j], a, x[j]);
            QD::mul(iy[j], b, ix[j]);
            QD::add(x[i], ix[j], y[j]);
            QD::sub(x[j], ix[i], iy[j]);

            QD::add(x[i_], x[j_], y[i_]);
            QD::add(ix[i_], ix[j_], iy[i_]);

            QD::sub(x[i_], x[j_], y[j_]);
            QD::sub(ix[i_], ix[j_], iy[j_]);

            QD::minus(a);
            QD::mul(y[j_], a, x[i_]);
            QD::mul(y[j_], b, ix[i_]);
            QD::mul(iy[j_], a, x[j_]);
            QD::mul(iy[j_], b, ix[j_]);
            QD::sub(x[i_], ix[j_], y[j_]);
            QD::add(x[j_], ix[i_], iy[j_]);

            i++;
            j++;
            i_--;
            j_--;
        }

        // y[0] = x[0] + x[n / 2]
        QD::add(x[i], x[j], y[i]);
        QD::add(ix[i], ix[j], iy[i]);
        // y[n / 2] = x[0] - x[n / 2]
        QD::sub(ix[i], ix[j], y[j]);
        QD::sub(x[j], x[i], iy[j]);
    }

    decimal_in_frequency(initN, harf_n, y, iy, x, ix, cos_table);
    decimal_in_frequency(initN, harf_n, y + harf_n, iy + harf_n, x + harf_n,
                         ix + harf_n, cos_table);

    for (i = 0; i < harf_n; i++) {
        QD::copy(y[i], x[2 * i]);
        QD::copy(iy[i], ix[2 * i]);
        QD::copy(y[i + harf_n], x[2 * i + 1]);
        QD::copy(iy[i + harf_n], ix[2 * i + 1]);
    }
}

}  // namespace FFT