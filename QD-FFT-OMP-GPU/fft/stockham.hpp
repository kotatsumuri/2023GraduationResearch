#pragma once
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace Stockham {
void fft(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd *y[], qd *iy[], qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)w[j * n / (2 * l)];
            double *b = (double *)iw[j * n / (2 * l)];
            for (uint64_t k = 0; k < m; k++) {
                butterfly((*x)[k + j * m], (*ix)[k + j * m],
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

void ifft(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd *y[], qd *iy[], qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)w[j * n / (2 * l)];
            double *b = (double *)iw[j * n / (2 * l)];
            for (uint64_t k = 0; k < m; k++) {
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
}  // namespace Stockham

void stockham(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    Stockham::fft(n, p, x, ix, &y, &iy, w, iw);
}

void inv_stockham(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    Stockham::ifft(n, p, x, ix, &y, &iy, w, iw);
    for (uint64_t i = 0; i < n; i++) {
        div_pwr2((*x)[i], n, (*x)[i]);
        div_pwr2((*ix)[i], n, (*ix)[i]);
    }
}