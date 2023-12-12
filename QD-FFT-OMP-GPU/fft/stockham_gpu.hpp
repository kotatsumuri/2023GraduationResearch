
#pragma once
#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace StockhamGPU {
inline void fft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }
    }
}

inline void fft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd y[], qd iy[], qd w[], qd iw[]) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            copy(x[i], y[i]);
            copy(ix[i], iy[i]);
        }
    }
}

inline void fft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[], Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    h2d_timer.start();
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        h2d_timer.stop();
        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            kernel_timer.stop();
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }
        d2h_timer.start();
    }
    d2h_timer.stop();
}

inline void fft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd y[], qd iy[], qd w[], qd iw[], Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    h2d_timer.start();
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        h2d_timer.stop();
        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            kernel_timer.stop();
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            copy(x[i], y[i]);
            copy(ix[i], iy[i]);
        }
        d2h_timer.start();
    }
    d2h_timer.stop();
}

inline void ifft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    inv_butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x[i], n, x[i]);
            div_pwr2(ix[i], n, ix[i]);
        }
    }
}

inline void ifft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], qd iw[]) {
    uint64_t n2 = n >> 1;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 2], iw[ : n / 2]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    double *a = (double *)w[j << (p - (1 + lp))];
                    double *b = (double *)iw[j << (p - (1 + lp))];
                    qd x0, ix0, x1, ix1, y0, iy0, y1, iy1;
                    copy(x[k + (j << mp)], x0);
                    copy(ix[k + (j << mp)], ix0);
                    copy(x[k + (j << mp) + n2], x1);
                    copy(ix[k + (j << mp) + n2], ix1);
                    inv_butterfly(x0, ix0, x1, ix1, y0, iy0, y1, iy1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(iy0, iy[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                    copy(iy1, iy[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            swap(&ix, &iy);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
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

void stockham_gpu(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], qd iw[], Timer &timer, Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    if (p & 1)
        StockhamGPU::fft_odd(n, p, *x, *ix, y, iy, w, iw, h2d_timer, d2h_timer, kernel_timer);
    else
        StockhamGPU::fft_even(n, p, *x, *ix, y, iy, w, iw, h2d_timer, d2h_timer, kernel_timer);
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