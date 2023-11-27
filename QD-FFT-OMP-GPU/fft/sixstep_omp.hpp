#pragma once
#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace SixStepOMP {
inline void fft(uint64_t n, uint64_t p, uint64_t u, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)w[j * n / (2 * l) * u];
            double *b = (double *)iw[j * n / (2 * l) * u];
            for (uint64_t k = 0; k < m; k++) {
                butterfly(x[k + j * m], ix[k + j * m],
                          x[k + j * m + l * m], ix[k + j * m + l * m],
                          y[k + 2 * j * m], iy[k + 2 * j * m],
                          y[k + 2 * j * m + m], iy[k + 2 * j * m + m],
                          a, b);
            }
        }
        swap(&x, &y);
        swap(&ix, &iy);
        l >>= 1;
        m <<= 1;
    }
}

inline void ifft(uint64_t n, uint64_t p, uint64_t u, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)w[j * n / (2 * l) * u];
            double *b = (double *)iw[j * n / (2 * l) * u];
            for (uint64_t k = 0; k < m; k++) {
                inv_butterfly(x[k + j * m], ix[k + j * m],
                              x[k + j * m + l * m], ix[k + j * m + l * m],
                              y[k + 2 * j * m], iy[k + 2 * j * m],
                              y[k + 2 * j * m + m], iy[k + 2 * j * m + m],
                              a, b);
            }
        }
        swap(&x, &y);
        swap(&ix, &iy);
        l >>= 1;
        m <<= 1;
    }
}

inline void twist(qd x, qd ix, qd y, qd iy, qd a, qd b) {
    qd tmp0, tmp1;
    mul(y, a, x);
    mul(y, b, tmp0);
    mul(iy, a, ix);
    mul(iy, b, tmp1);
    add(x, tmp1, x);
    sub(ix, tmp0, ix);
}

inline void inv_twist(qd x, qd ix, qd y, qd iy, qd a, qd b) {
    qd tmp0, tmp1;
    mul(y, a, x);
    mul(y, b, tmp0);
    mul(iy, a, ix);
    mul(iy, b, tmp1);
    sub(x, tmp1, x);
    add(ix, tmp0, ix);
}

void sixstep(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd *y[], qd *iy[], qd w[], qd iw[]) {
    uint64_t p1 = p >> 1;
    uint64_t p2 = p - p1;
    uint64_t n1 = 1ull << p1;
    uint64_t n2 = 1ull << p2;
    uint64_t u1 = n >> p1;
    uint64_t u2 = n >> p2;

#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < n1; i++) {
        for (uint64_t j = 0; j < n2; j++) {
            copy((*x)[i * n2 + j], (*y)[j * n1 + i]);
            copy((*ix)[i * n2 + j], (*iy)[j * n1 + i]);
        }
    }

#pragma omp parallel for
    for (uint64_t j = 0; j < n2; j++) {
        fft(n1, p1, u1, (*y) + j * n1, (*iy) + j * n1, (*x) + j * n1, (*ix) + j * n1, w, iw);
    }

    if (p1 & 1) {
        swap(x, y);
        swap(ix, iy);
    }

#pragma omp parallel for collapse(2)
    for (uint64_t j = 0; j < n2; j++) {
        for (uint64_t i = 0; i < n1; i++) {
            double *a = (double *)w[i * j];
            double *b = (double *)iw[i * j];
            twist((*x)[i * n2 + j], (*ix)[i * n2 + j], (*y)[j * n1 + i], (*iy)[j * n1 + i], a, b);
        }
    }

#pragma omp parallel for
    for (uint64_t i = 0; i < n1; i++) {
        fft(n2, p2, u2, (*x) + i * n2, (*ix) + i * n2, (*y) + i * n2, (*iy) + i * n2, w, iw);
    }

    if (p2 & 1) {
        swap(x, y);
        swap(ix, iy);
    }

#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < n1; i++) {
        for (uint64_t j = 0; j < n2; j++) {
            copy((*x)[i * n2 + j], (*y)[j * n1 + i]);
            copy((*ix)[i * n2 + j], (*iy)[j * n1 + i]);
        }
    }
}

void inv_sixstep(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd *y[], qd *iy[], qd w[], qd iw[]) {
    uint64_t p1 = p >> 1;
    uint64_t p2 = p - p1;
    uint64_t n1 = 1ull << p1;
    uint64_t n2 = 1ull << p2;
    uint64_t u1 = n >> p1;
    uint64_t u2 = n >> p2;

#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < n1; i++) {
        for (uint64_t j = 0; j < n2; j++) {
            copy((*x)[i * n2 + j], (*y)[j * n1 + i]);
            copy((*ix)[i * n2 + j], (*iy)[j * n1 + i]);
        }
    }

#pragma omp parallel for
    for (uint64_t j = 0; j < n2; j++) {
        ifft(n1, p1, u1, (*y) + j * n1, (*iy) + j * n1, (*x) + j * n1, (*ix) + j * n1, w, iw);
    }

    if (p1 & 1) {
        swap(x, y);
        swap(ix, iy);
    }

#pragma omp parallel for collapse(2)
    for (uint64_t j = 0; j < n2; j++) {
        for (uint64_t i = 0; i < n1; i++) {
            double *a = (double *)w[i * j];
            double *b = (double *)iw[i * j];
            inv_twist((*x)[i * n2 + j], (*ix)[i * n2 + j], (*y)[j * n1 + i], (*iy)[j * n1 + i], a, b);
        }
    }

#pragma omp parallel for
    for (uint64_t i = 0; i < n1; i++) {
        ifft(n2, p2, u2, (*x) + i * n2, (*ix) + i * n2, (*y) + i * n2, (*iy) + i * n2, w, iw);
    }

    if (p2 & 1) {
        swap(x, y);
        swap(ix, iy);
    }

#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < n1; i++) {
        for (uint64_t j = 0; j < n2; j++) {
            div_pwr2((*x)[i * n2 + j], n, (*y)[j * n1 + i]);
            div_pwr2((*ix)[i * n2 + j], n, (*iy)[j * n1 + i]);
        }
    }
}
}  // namespace SixStepOMP

void sixstep_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    SixStepOMP::sixstep(n, p, x, ix, &y, &iy, w, iw);
    swap(x, &y);
    swap(ix, &iy);
    free(y);
    free(iy);
}

void sixstep_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[], Timer &timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    SixStepOMP::sixstep(n, p, x, ix, &y, &iy, w, iw);
    timer.stop();
    swap(x, &y);
    swap(ix, &iy);
    free(y);
    free(iy);
}

void inv_sixstep_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    SixStepOMP::inv_sixstep(n, p, x, ix, &y, &iy, w, iw);
    swap(x, &y);
    swap(ix, &iy);
    free(y);
    free(iy);
}

void inv_sixstep_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[], Timer &timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    SixStepOMP::inv_sixstep(n, p, x, ix, &y, &iy, w, iw);
    timer.stop();
    swap(x, &y);
    swap(ix, &iy);
    free(y);
    free(iy);
}