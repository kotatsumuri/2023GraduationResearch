#pragma once
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace DIF {
void fft(uint64_t n, uint64_t m, uint64_t u, double *x[], double *ix[], double *y[], double *iy[], qd w[], qd iw[]) {
    if (m <= 1)
        return;

    uint64_t m2 = m >> 1;

    for (uint64_t i = 0, j = m2, k = 0; i < m2; i++, j++, k += u) {
        double *a = (double *)w[k];
        double *b = (double *)iw[k];
        butterfly(x[i], ix[i], x[j], ix[j], y[i], iy[i], y[j], iy[j], a, b);
    }

    uint64_t u2 = u << 1;
    fft(n, m2, u2, y, iy, x, ix, w, iw);
    fft(n, m2, u2, y + m2, iy + m2, x + m2, ix + m2, w, iw);

    for (uint64_t i = 0, j = 0, k = m2, l = 1; i < m2; i++, j += 2, k++, l += 2) {
        swap(y + i, x + j);
        swap(iy + i, ix + j);
        swap(y + k, x + l);
        swap(iy + k, ix + l);
    }
}

void ifft(uint64_t n, uint64_t m, uint64_t u, double *x[], double *ix[], double *y[], double *iy[], qd w[], qd iw[]) {
    if (m <= 1)
        return;

    uint64_t m2 = m >> 1;

    for (uint64_t i = 0, j = m2, k = 0; i < m2; i++, j++, k += u) {
        double *a = (double *)w[k];
        double *b = (double *)iw[k];
        inv_butterfly(x[i], ix[i], x[j], ix[j], y[i], iy[i], y[j], iy[j], a, b);
    }

    uint64_t u2 = u << 1;
    ifft(n, m2, u2, y, iy, x, ix, w, iw);
    ifft(n, m2, u2, y + m2, iy + m2, x + m2, ix + m2, w, iw);

    for (uint64_t i = 0, j = 0, k = m2, l = 1; i < m2; i++, j += 2, k++, l += 2) {
        swap(y + i, x + j);
        swap(iy + i, ix + j);
        swap(y + k, x + l);
        swap(iy + k, ix + l);
    }
}
}  // namespace DIF

void dif(uint64_t n, qd x[], qd ix[], qd w[], qd iw[]) {
    qd y[n];
    qd iy[n];

    double *x_[n];
    double *ix_[n];
    double *y_[n];
    double *iy_[n];

    for (uint64_t i = 0; i < n; i++) {
        x_[i]  = x[i];
        ix_[i] = ix[i];
        y_[i]  = y[i];
        iy_[i] = iy[i];
    }

    DIF::fft(n, n, 1, x_, ix_, y_, iy_, w, iw);

    qd z[n];
    qd iz[n];

    copy_vector(n, x_, z);
    copy_vector(n, ix_, iz);

    copy_vector(n, z, x);
    copy_vector(n, iz, ix);
}

void inv_dif(uint64_t n, qd x[], qd ix[], qd w[], qd iw[]) {
    qd y[n];
    qd iy[n];

    double *x_[n];
    double *ix_[n];
    double *y_[n];
    double *iy_[n];

    for (uint64_t i = 0; i < n; i++) {
        x_[i]  = x[i];
        ix_[i] = ix[i];
        y_[i]  = y[i];
        iy_[i] = iy[i];
    }

    DIF::ifft(n, n, 1, x_, ix_, y_, iy_, w, iw);

    qd z[n];
    qd iz[n];

    copy_vector(n, x_, z);
    copy_vector(n, ix_, iz);

    for (uint64_t i = 0; i < n; i++) {
        div_pwr2(z[i], n, x[i]);
        div_pwr2(iz[i], n, ix[i]);
    }
}