#pragma once
#include <iostream>

#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace StockhamGPUN4 {
inline void fft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, const qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    uint64_t idx_n4 = j >> (lp - 1);
                    uint64_t idx    = j << (p - (1 + lp));
                    uint64_t idx_a  = idx & (n4 - 1);
                    if (idx_n4)
                        idx_a = n4 - idx_a;
                    double *a = (double *)w[idx_a];
                    qd a_;
                    if (idx_n4) {
                        minus(a, a_);
                        a = (double *)a_;
                    }
                    idx_a     = n4 - idx_a;
                    double *b = (double *)w[idx_a];
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

inline void fft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, const qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    uint64_t idx_n4 = j >> (lp - 1);
                    uint64_t idx    = j << (p - (1 + lp));
                    uint64_t idx_a  = idx & (n4 - 1);
                    if (idx_n4)
                        idx_a = n4 - idx_a;
                    double *a = (double *)w[idx_a];
                    qd a_;
                    if (idx_n4) {
                        minus(a, a_);
                        a = (double *)a_;
                    }
                    idx_a     = n4 - idx_a;
                    double *b = (double *)w[idx_a];
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

inline void fft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[], Timer &h2dtimer, Timer &d2htimer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    qd *w_n2    = (qd *)calloc(n2 + 1, sizeof(qd));
    qd *iw_n2   = (qd *)calloc(n2 + 1, sizeof(qd));
    h2dtimer.start();
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n], w_n2[ : n2 + 1], iw_n2[ : n2 + 1])
    {
        h2dtimer.stop();

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n4 + 1; i++) {
            copy(w[i], w_n2[i]);
            minus(w[i], w_n2[n2 - i]);
            copy(w[i], iw_n2[n4 - i]);
            copy(w[i], iw_n2[n4 + i]);
        }

        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    // uint64_t idx_n4 = j >> (lp - 1);
                    // uint64_t idx    = (j << (p - (1 + lp))) & (n4 - 1);
                    // if (idx_n4)
                    //     idx = n4 - idx;
                    // double *a = (double *)w[idx];
                    // qd a_;
                    // if (idx_n4) {
                    //     minus(a, a_);
                    //     a = (double *)a_;
                    // }
                    // double *b = (double *)w[n4 - idx];
                    double *a = (double *)w_n2[j << (p - (1 + lp))];
                    double *b = (double *)iw_n2[j << (p - (1 + lp))];
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
        d2htimer.start();
    }
    d2htimer.stop();
}

inline void fft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd y[], qd iy[], qd w[], Timer &h2dtimer, Timer &d2htimer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    qd *w_n2    = (qd *)calloc(n2 + 1, sizeof(qd));
    qd *iw_n2   = (qd *)calloc(n2 + 1, sizeof(qd));
    h2dtimer.start();
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n], w_n2[ : n2 + 1], iw_n2[ : n2 + 1])
    {
        h2dtimer.stop();

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n4 + 1; i++) {
            copy(w[i], w_n2[i]);
            minus(w[i], w_n2[n2 - i]);
            copy(w[i], iw_n2[n4 - i]);
            copy(w[i], iw_n2[n4 + i]);
        }

        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    // uint64_t idx_n4 = j >> (lp - 1);
                    // uint64_t idx    = (j << (p - (1 + lp))) & (n4 - 1);
                    // if (idx_n4)
                    //     idx = n4 - idx;
                    // double *a = (double *)w[idx];
                    // qd a_;
                    // if (idx_n4) {
                    //     minus(a, a_);
                    //     a = (double *)a_;
                    // }
                    // double *b = (double *)w[n4 - idx];
                    double *a = (double *)w_n2[j << (p - (1 + lp))];
                    double *b = (double *)iw_n2[j << (p - (1 + lp))];
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
        d2htimer.start();
    }
    d2htimer.stop();
}

inline void ifft_even(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    uint64_t idx_n4 = j >> (lp - 1);
                    uint64_t idx    = j << (p - (1 + lp));
                    uint64_t idx_a  = idx & (n4 - 1);
                    if (idx_n4)
                        idx_a = n4 - idx_a;
                    double *a = (double *)w[idx_a];
                    qd a_;
                    if (idx_n4) {
                        minus(a, a_);
                        a = (double *)a_;
                    }
                    idx_a     = n4 - idx_a;
                    double *b = (double *)w[idx_a];
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

inline void ifft_odd(uint64_t n, uint64_t p, qd *x, qd *ix, qd *y, qd *iy, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n], ix[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n], iy[ : n])
    {
        for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < l; j++) {
                for (uint64_t k = 0; k < m; k++) {
                    uint64_t idx_n4 = j >> (lp - 1);
                    uint64_t idx    = j << (p - (1 + lp));
                    uint64_t idx_a  = idx & (n4 - 1);
                    if (idx_n4)
                        idx_a = n4 - idx_a;
                    double *a = (double *)w[idx_a];
                    qd a_;
                    if (idx_n4) {
                        minus(a, a_);
                        a = (double *)a_;
                    }
                    idx_a     = n4 - idx_a;
                    double *b = (double *)w[idx_a];
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
}  // namespace StockhamGPUN4

void stockham_gpu_n4(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    if (p & 1)
        StockhamGPUN4::fft_odd(n, p, *x, *ix, y, iy, w);
    else
        StockhamGPUN4::fft_even(n, p, *x, *ix, y, iy, w);
    free(y);
    free(iy);
}

void stockham_gpu_n4(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[], Timer &timer, Timer &h2dtimer, Timer &d2htimer, Timer &kernel_timer) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    timer.start();
    if (p & 1)
        StockhamGPUN4::fft_odd(n, p, *x, *ix, y, iy, w, h2dtimer, d2htimer, kernel_timer);
    else
        StockhamGPUN4::fft_even(n, p, *x, *ix, y, iy, w, h2dtimer, d2htimer, kernel_timer);
    timer.stop();
    free(y);
    free(iy);
}

void inv_stockham_gpu_n4(uint64_t n, uint64_t p, qd *x[], qd *ix[], qd w[]) {
    qd *y  = (qd *)calloc(n, sizeof(qd));
    qd *iy = (qd *)calloc(n, sizeof(qd));
    if (p & 1)
        StockhamGPUN4::ifft_odd(n, p, *x, *ix, y, iy, w);
    else
        StockhamGPUN4::ifft_even(n, p, *x, *ix, y, iy, w);
    free(y);
    free(iy);
}