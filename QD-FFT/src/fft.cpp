#include "fft.hpp"

#include <cmath>
#include <iostream>
#include <queue>

#include "qd.hpp"

namespace FFT {
void make_cos_table(long long int N, QD::qd cos_table[]) {
    QD::init(cos_table[0], 1.0);
    if (N <= 2)
        return;
    QD::init(cos_table[N / 4], 0.0);
    if (N <= 4)
        return;
    QD::sqrt(0.5, cos_table[N / 8]);
    if (N <= 8)
        return;
    std::queue<int> que;
    QD::qd h;
    long long int i, j, k;
    que.push(N / 8);
    while (!que.empty()) {
        i = que.front();
        j = i / 2;
        k = N / 4 - j;
        que.pop();
        QD::mul_pwr2(cos_table[i], 2, h);
        QD::add(h, 2.0, cos_table[k]);
        QD::sqrt(cos_table[k], cos_table[j]);
        QD::mul_pwr2(cos_table[j], 0.5, cos_table[j]);
        QD::sub(2.0, h, cos_table[k]);
        QD::sqrt(cos_table[k], h);
        QD::mul_pwr2(h, 0.5, cos_table[k]);
        if (!(j & 1)) {
            que.push(j);
            que.push(k);
        }
    }
}

void make_cos_table(long long int N, double cos_table[]) {
    cos_table[0] = 1.0;
    if (N <= 2)
        return;
    cos_table[N / 4] = 0.0;
    if (N <= 4)
        return;
    cos_table[N / 8] = std::sqrt(0.5);
    if (N <= 8)
        return;
    std::queue<int> que;
    long long int i, j, k;
    double h;
    que.push(N / 8);
    while (!que.empty()) {
        i = que.front();
        j = i / 2;
        k = N / 4 - j;
        que.pop();
        h            = cos_table[i] * 2;
        cos_table[j] = 0.5 * std::sqrt(2.0 + h);
        cos_table[k] = 0.5 * std::sqrt(2.0 - h);
        if (!(j & 1)) {
            que.push(j);
            que.push(k);
        }
    }
}

void decimation_in_frequency(long long int initN, long long int n, QD::qd x[],
                             QD::qd ix[], QD::qd y[], QD::qd iy[],
                             QD::qd cos_table[]) {
    const long long int harf_n = n / 2;
    // y[0] = x[0] + x[n / 2]
    QD::add(x[0], x[harf_n], y[0]);
    QD::add(ix[0], ix[harf_n], iy[0]);
    // y[n / 2] = x[0] - x[n / 2]
    QD::sub(x[0], x[harf_n], y[harf_n]);
    QD::sub(ix[0], ix[harf_n], iy[harf_n]);

    long long int i = 1;
    if (n > 2) {
        long long int j                 = 1 + harf_n;
        long long int i_                = harf_n - 1;
        long long int j_                = n - 1;
        const long long int harf_harf_n = harf_n / 2;
        const long long int initN_n     = initN / n;
        const long long int initN_4     = initN / 4;

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

            QD::sub(x[j_], x[i_], y[j_]);
            QD::sub(ix[j_], ix[i_], iy[j_]);

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

        // y[n / 4] = x[n / 4] + x[3 * n / 4]
        QD::add(x[i], x[j], y[i]);
        QD::add(ix[i], ix[j], iy[i]);
        // y[3 * n / 4] = x[n / 4] - x[3 * n / 4]
        QD::sub(ix[i], ix[j], y[j]);
        QD::sub(x[j], x[i], iy[j]);

        decimation_in_frequency(initN, harf_n, y, iy, x, ix, cos_table);
        decimation_in_frequency(initN, harf_n, y + harf_n, iy + harf_n,
                                x + harf_n, ix + harf_n, cos_table);

        for (i = 0; i < harf_n; i++) {
            std::swap(y[i], x[2 * i]);
            std::swap(iy[i], ix[2 * i]);
            std::swap(y[i + harf_n], x[2 * i + 1]);
            std::swap(iy[i + harf_n], ix[2 * i + 1]);
        }
    } else {
        std::swap(y[0], x[0]);
        std::swap(iy[0], ix[0]);
        std::swap(y[1], x[1]);
        std::swap(iy[1], ix[1]);
    }
}

void decimation_in_frequency(long long int initN, long long int n, double x[],
                             double ix[], double y[], double iy[],
                             double cos_table[]) {
    const long long int harf_n = n / 2;
    // y[0] = x[0] + x[n / 2]
    y[0]  = x[0] + x[harf_n];
    iy[0] = ix[0] + ix[harf_n];
    // y[n / 2] = x[0] - x[n / 2]
    y[harf_n]  = x[0] - x[harf_n];
    iy[harf_n] = ix[0] - ix[harf_n];

    long long int i = 1;
    if (n > 2) {
        long long int j                 = 1 + harf_n;
        long long int i_                = harf_n - 1;
        long long int j_                = n - 1;
        const long long int harf_harf_n = harf_n / 2;
        const long long int initN_n     = initN / n;
        const long long int initN_4     = initN / 4;

        for (; i < harf_harf_n;) {
            y[i]  = x[i] + x[j];
            iy[i] = ix[i] + ix[j];

            y[j]  = x[i] - x[j];
            iy[j] = ix[i] - ix[j];

            double a = cos_table[i * initN_n];
            double b = cos_table[initN_4 - i * initN_n];

            x[i]  = y[j] * a;
            ix[i] = y[j] * b;
            x[j]  = iy[j] * a;
            ix[j] = iy[j] * b;
            y[j]  = x[i] + ix[j];
            iy[j] = x[j] - ix[i];

            y[i_]  = x[i_] + x[j_];
            iy[i_] = ix[i_] + ix[j_];

            y[j_]  = x[j_] - x[i_];
            iy[j_] = ix[j_] - ix[i_];

            x[i_]  = y[j_] * a;
            ix[i_] = y[j_] * b;
            x[j_]  = iy[j_] * a;
            ix[j_] = iy[j_] * b;
            y[j_]  = x[i_] - ix[j_];
            iy[j_] = x[j_] + ix[i_];

            i++;
            j++;
            i_--;
            j_--;
        }

        // y[n / 4] = x[n / 4] + x[3 * n / 4]
        y[i]  = x[i] + x[j];
        iy[i] = ix[i] + ix[j];
        y[j]  = ix[i] - ix[j];
        iy[j] = x[j] - x[i];

        decimation_in_frequency(initN, harf_n, y, iy, x, ix, cos_table);
        decimation_in_frequency(initN, harf_n, y + harf_n, iy + harf_n,
                                x + harf_n, ix + harf_n, cos_table);

        for (i = 0; i < harf_n; i++) {
            x[2 * i]      = y[i];
            ix[2 * i]     = iy[i];
            x[2 * i + 1]  = y[i + harf_n];
            ix[2 * i + 1] = iy[i + harf_n];
        }
    } else {
        x[0]  = y[0];
        ix[0] = iy[0];
        x[1]  = y[1];
        ix[1] = iy[1];
    }
}

void inv_decimation_in_frequency(long long int initN, long long int n,
                                 double x[], double ix[], double y[],
                                 double iy[], double cos_table[]) {
    const long long int harf_n = n / 2;
    // y[0] = x[0] + x[n / 2]
    y[0]  = x[0] + x[harf_n];
    iy[0] = ix[0] + ix[harf_n];
    // y[n / 2] = x[0] - x[n / 2]
    y[harf_n]  = x[0] - x[harf_n];
    iy[harf_n] = ix[0] - ix[harf_n];

    long long int i = 1;
    if (n > 2) {
        long long int j                 = 1 + harf_n;
        long long int i_                = harf_n - 1;
        long long int j_                = n - 1;
        const long long int harf_harf_n = harf_n / 2;
        const long long int initN_n     = initN / n;
        const long long int initN_4     = initN / 4;

        for (; i < harf_harf_n;) {
            y[i]  = x[i] + x[j];
            iy[i] = ix[i] + ix[j];

            y[j]  = x[i] - x[j];
            iy[j] = ix[i] - ix[j];

            double a = cos_table[i * initN_n];
            double b = cos_table[initN_4 - i * initN_n];

            x[i]  = y[j] * a;
            ix[i] = y[j] * b;
            x[j]  = iy[j] * a;
            ix[j] = iy[j] * b;
            y[j]  = x[i] - ix[j];
            iy[j] = x[j] + ix[i];

            y[i_]  = x[i_] + x[j_];
            iy[i_] = ix[i_] + ix[j_];

            y[j_]  = x[j_] - x[i_];
            iy[j_] = ix[j_] - ix[i_];

            x[i_]  = y[j_] * a;
            ix[i_] = y[j_] * b;
            x[j_]  = iy[j_] * a;
            ix[j_] = iy[j_] * b;
            y[j_]  = x[i_] + ix[j_];
            iy[j_] = x[j_] - ix[i_];

            i++;
            j++;
            i_--;
            j_--;
        }

        // y[n / 4] = x[n / 4] + x[3 * n / 4]
        y[i]  = x[i] + x[j];
        iy[i] = ix[i] + ix[j];
        y[j]  = ix[j] - ix[i];
        iy[j] = x[i] - x[j];

        inv_decimation_in_frequency(initN, harf_n, y, iy, x, ix, cos_table);
        inv_decimation_in_frequency(initN, harf_n, y + harf_n, iy + harf_n,
                                    x + harf_n, ix + harf_n, cos_table);

        for (i = 0; i < harf_n; i++) {
            std::swap(y[i], x[2 * i]);
            std::swap(iy[i], ix[2 * i]);
            std::swap(y[i + harf_n], x[2 * i + 1]);
            std::swap(iy[i + harf_n], ix[2 * i + 1]);
        }
    } else {
        x[0]  = y[0] / initN;
        ix[0] = iy[0] / initN;
        x[1]  = y[1] / initN;
        ix[1] = iy[1] / initN;
    }
}

void inv_decimation_in_frequency(long long int initN, long long int n,
                                 QD::qd x[], QD::qd ix[], QD::qd y[],
                                 QD::qd iy[], QD::qd cos_table[]) {
    const long long int harf_n = n / 2;
    // y[0] = x[0] + x[n / 2]
    QD::add(x[0], x[harf_n], y[0]);
    QD::add(ix[0], ix[harf_n], iy[0]);
    // y[n / 2] = x[0] - x[n / 2]
    QD::sub(x[0], x[harf_n], y[harf_n]);
    QD::sub(ix[0], ix[harf_n], iy[harf_n]);

    long long int i = 1;
    if (n > 2) {
        long long int j                 = 1 + harf_n;
        long long int i_                = harf_n - 1;
        long long int j_                = n - 1;
        const long long int harf_harf_n = harf_n / 2;
        const long long int initN_n     = initN / n;
        const long long int initN_4     = initN / 4;

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
            QD::sub(x[i], ix[j], y[j]);
            QD::add(x[j], ix[i], iy[j]);

            QD::add(x[i_], x[j_], y[i_]);
            QD::add(ix[i_], ix[j_], iy[i_]);

            QD::sub(x[j_], x[i_], y[j_]);
            QD::sub(ix[j_], ix[i_], iy[j_]);

            QD::mul(y[j_], a, x[i_]);
            QD::mul(y[j_], b, ix[i_]);
            QD::mul(iy[j_], a, x[j_]);
            QD::mul(iy[j_], b, ix[j_]);
            QD::add(x[i_], ix[j_], y[j_]);
            QD::sub(x[j_], ix[i_], iy[j_]);

            i++;
            j++;
            i_--;
            j_--;
        }

        // y[n / 4] = x[n / 4] + x[3 * n / 4]
        QD::add(x[i], x[j], y[i]);
        QD::add(ix[i], ix[j], iy[i]);
        // y[3 * n / 4] = x[n / 4] - x[3 * n / 4]
        QD::sub(ix[j], ix[i], y[j]);
        QD::sub(x[i], x[j], iy[j]);

        inv_decimation_in_frequency(initN, harf_n, y, iy, x, ix, cos_table);
        inv_decimation_in_frequency(initN, harf_n, y + harf_n, iy + harf_n,
                                    x + harf_n, ix + harf_n, cos_table);

        for (i = 0; i < harf_n; i++) {
            std::swap(y[i], x[2 * i]);
            std::swap(iy[i], ix[2 * i]);
            std::swap(y[i + harf_n], x[2 * i + 1]);
            std::swap(iy[i + harf_n], ix[2 * i + 1]);
        }
    } else {
        double inv_initN = 1.0 / initN;
        QD::mul_pwr2(y[0], inv_initN, x[0]);
        QD::mul_pwr2(iy[0], inv_initN, ix[0]);
        QD::mul_pwr2(y[1], inv_initN, x[1]);
        QD::mul_pwr2(iy[1], inv_initN, ix[1]);
    }
}

void DFT(int N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
         QD::qd cos_table[]) {
    if (N == 2) {
        QD::add(x[0], x[1], y[0]);
        QD::add(ix[0], ix[1], iy[0]);
        QD::sub(x[0], x[1], y[0]);
        QD::sub(ix[0], ix[1], iy[0]);
        return;
    }
    QD::qd _cos_table[N];
    QD::qd _sin_table[N];
    QD::init(_cos_table[0], 1);
    QD::zero(_cos_table[N / 4]);
    QD::init(_cos_table[N / 2], -1);
    QD::zero(_cos_table[3 * N / 4]);
    QD::zero(_sin_table[0]);
    QD::init(_sin_table[N / 4], 1);
    QD::zero(_sin_table[N / 2]);
    QD::init(_sin_table[3 * N / 4], -1);
    for (int i = 1; i < N / 4; i++) {
        QD::copy(cos_table[i], _cos_table[i]);
        QD::minus(cos_table[N / 4 - i], _cos_table[i + N / 4]);
        QD::copy(_cos_table[i], _cos_table[N - i]);
        QD::copy(_cos_table[i + N / 4], _cos_table[3 * N / 4 - i]);
        QD::copy(cos_table[N - i], _sin_table[i]);
        QD::copy(_sin_table[i], _sin_table[N / 2 - i]);
        QD::minus(_sin_table[i], _sin_table[N / 2 + i]);
        QD::copy(_sin_table[N / 2 + i], _sin_table[N - i]);
    }

    for (int i = 0; i < N; i++) {
        QD::zero(y[i]);
    }

    QD::qd t;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            QD::mul(x[j], _cos_table[i * j % N], t);
            QD::add(y[i], t, y[i]);
            QD::mul(ix[j], _sin_table[i * j % N], t);
            QD::add(y[i], t, y[i]);
            QD::mul(x[j], _sin_table[i * j % N], t);
            QD::sub(iy[i], t, iy[i]);
            QD::mul(ix[j], _cos_table[i * j % N], t);
            QD::add(iy[i], t, iy[i]);
        }
    }

    for (int i = 0; i < N; i++) {
        std::swap(x[i], y[i]);
    }
}

void inv_DFT(int N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
             QD::qd cos_table[]) {
    if (N == 2) {
        QD::add(x[0], x[1], y[0]);
        QD::add(ix[0], ix[1], iy[0]);
        QD::sub(x[0], x[1], y[1]);
        QD::sub(ix[0], ix[1], iy[1]);
        return;
    }
    QD::qd _cos_table[N];
    QD::qd _sin_table[N];
    QD::init(_cos_table[0], 1);
    QD::zero(_cos_table[N / 4]);
    QD::init(_cos_table[N / 2], -1);
    QD::zero(_cos_table[3 * N / 4]);
    QD::zero(_sin_table[0]);
    QD::init(_sin_table[N / 4], 1);
    QD::zero(_sin_table[N / 2]);
    QD::init(_sin_table[3 * N / 4], -1);
    for (int i = 1; i < N / 4; i++) {
        QD::copy(cos_table[i], _cos_table[i]);
        QD::minus(cos_table[N / 4 - i], _cos_table[i + N / 4]);
        QD::copy(_cos_table[i], _cos_table[N - i]);
        QD::copy(_cos_table[i + N / 4], _cos_table[3 * N / 4 - i]);
        QD::copy(cos_table[N - i], _sin_table[i]);
        QD::copy(_sin_table[i], _sin_table[N / 2 - i]);
        QD::minus(_sin_table[i], _sin_table[N / 2 + i]);
        QD::copy(_sin_table[N / 2 + i], _sin_table[N - i]);
    }

    for (int i = 0; i < N; i++) {
        QD::zero(y[i]);
    }

    QD::qd t;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            QD::mul(x[j], _cos_table[i * j % N], t);
            QD::add(y[i], t, y[i]);
            QD::mul(ix[j], _sin_table[i * j % N], t);
            QD::sub(y[i], t, y[i]);
            QD::mul(x[j], _sin_table[i * j % N], t);
            QD::add(iy[i], t, iy[i]);
            QD::mul(ix[j], _cos_table[i * j % N], t);
            QD::add(iy[i], t, iy[i]);
        }
        QD::mul_pwr2(y[i], 1.0 / N, y[i]);
    }

    for (int i = 0; i < N; i++) {
        std::swap(x[i], y[i]);
    }
}

}  // namespace FFT