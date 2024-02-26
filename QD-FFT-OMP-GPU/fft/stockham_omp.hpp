#pragma once
#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace StockhamOMP {
void fft(uint64_t n, uint64_t p, qd x[], qd ix[], qd y[], qd iy[], qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
#pragma omp parallel for collapse(2)
        for (uint64_t j = 0; j < l; j++) {
            for (uint64_t k = 0; k < m; k++) {
                double *a = (double *)w[j * n / (2 * l)];
                double *b = (double *)iw[j * n / (2 * l)];
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

void ifft(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd *y[], qd *iy[], qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
#pragma omp parallel for collapse(2)
        for (uint64_t j = 0; j < l; j++) {
            for (uint64_t k = 0; k < m; k++) {
                double *a = (double *)w[j * n / (2 * l)];
                double *b = (double *)iw[j * n / (2 * l)];
                inv_butterfly((*x)[k + j * m], (*ix)[k + j * m],
                              (*x)[k + j * m + l * m], (*ix)[k + j * m + l * m],
                              (*y)[k + 2 * j * m], (*iy)[k + 2 * j * m],
                              (*y)[k + 2 * j * m + m], (*iy)[k + 2 * j * m + m],
                              a, b);
            }
        }
        swap(x, y);
        swap(ix, iy);
        l >>= 1;
        m <<= 1;
    }
}
}  // namespace StockhamOMP

void stockham_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    StockhamOMP::fft(n, p, *x, *ix, y, iy, w, iw);
    if (p & 1) {
        swap(x, &y);
        swap(ix, &iy);
    }

    free(y);
    free(iy);
}

void stockham_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[], Timer &timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    StockhamOMP::fft(n, p, *x, *ix, y, iy, w, iw);
    if (p & 1) {
        swap(x, &y);
        swap(ix, &iy);
    }
    timer.stop();
    free(y);
    free(iy);
}

void inv_stockham_omp(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    StockhamOMP ::ifft(n, p, x, ix, &y, &iy, w, iw);
    for (uint64_t i = 0; i < n; i++) {
        div_pwr2((*x)[i], n, (*x)[i]);
        div_pwr2((*ix)[i], n, (*ix)[i]);
    }
    free(y);
    free(iy);
}
