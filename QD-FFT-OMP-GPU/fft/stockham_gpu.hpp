#pragma once
#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace StockhamGPU {
inline void fft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n], iw[ : n]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
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
}

inline void fft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd y[], qd iy[], qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n], iw[ : n]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
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

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            copy(x[i], y[i]);
            copy(ix[i], iy[i]);
        }
    }
}

inline void ifft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n], iw[ : n]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j * n / (2 * l)];
                    double *b = (double *)iw[j * n / (2 * l)];
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

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x[i], n, x[i]);
            div_pwr2(ix[i], n, ix[i]);
        }
    }
}

inline void ifft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t l = n >> 1;
    uint64_t m = 1;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n], iw[ : n]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j * n / (2 * l)];
                    double *b = (double *)iw[j * n / (2 * l)];
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

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x[i], n, y[i]);
            div_pwr2(ix[i], n, iy[i]);
        }
    }
}
}  // namespace StockhamGPU

void stockham_gpu(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    if (p & 1)
        StockhamGPU::fft_odd(n, p, *x, *ix, y, iy, w, iw);
    else
        StockhamGPU::fft_even(n, p, *x, *ix, y, iy, w, iw);
    free(y);
    free(iy);
}

void stockham_gpu(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[], Timer &timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    if (p & 1)
        StockhamGPU::fft_odd(n, p, *x, *ix, y, iy, w, iw);
    else
        StockhamGPU::fft_even(n, p, *x, *ix, y, iy, w, iw);
    timer.stop();
    free(y);
    free(iy);
}

void inv_stockham_gpu(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    if (p & 1)
        StockhamGPU::ifft_odd(n, p, *x, *ix, y, iy, w, iw);
    else
        StockhamGPU::ifft_even(n, p, *x, *ix, y, iy, w, iw);
    free(y);
    free(iy);
}