#pragma once
#include "../bench_util/Timer.hpp"
#include "../qd/qd.hpp"
#include "butterfly.hpp"
#include "fft_util.hpp"

namespace StockhamGPUN4Complex {
inline void fft_even(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
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
                    qd_complex x0, x1, y0, y1;
                    copy(x[k + (j << mp)], x0);
                    copy(x[k + (j << mp) + n2], x1);
                    butterfly(x0, x1, y0, y1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }
    }
}

inline void fft_odd(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
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
                    qd_complex x0, x1, y0, y1;
                    copy(x[k + (j << mp)], x0);
                    copy(x[k + (j << mp) + n2], x1);
                    butterfly(x0, x1, y0, y1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                }
            }
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            copy(x[i].re, y[i].re);
            copy(x[i].im, y[i].im);
        }
    }
}

inline void fft_even(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[], Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    h2d_timer.start();
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
    {
        h2d_timer.stop();
        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
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
                    qd_complex x0, x1, y0, y1;
                    copy(x[k + (j << mp)], x0);
                    copy(x[k + (j << mp) + n2], x1);
                    butterfly(x0, x1, y0, y1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                }
            }
            kernel_timer.stop();
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }
        d2h_timer.start();
    }
    d2h_timer.stop();
}

inline void fft_odd(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[], Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
    h2d_timer.start();
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
    {
        h2d_timer.stop();
        for (uint64_t t = 0; t < p; t++) {
            kernel_timer.start();
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
                    qd_complex x0, x1, y0, y1;
                    copy(x[k + (j << mp)], x0);
                    copy(x[k + (j << mp) + n2], x1);
                    butterfly(x0, x1, y0, y1, a, b);
                    copy(y0, y[k + (j << (mp + 1))]);
                    copy(y1, y[k + (j << (mp + 1)) + m]);
                }
            }
            kernel_timer.stop();
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            copy(x[i].re, y[i].re);
            copy(x[i].im, y[i].im);
        }
        d2h_timer.start();
    }
    d2h_timer.stop();
}

inline void ifft_even(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
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
                    inv_butterfly(x[k + (j << mp)],
                                  x[k + (j << mp) + n2],
                                  y[k + (j << (mp + 1))],
                                  y[k + (j << (mp + 1)) + m],
                                  a, b);
                }
            }
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x[i].re, n, x[i].re);
            div_pwr2(x[i].im, n, x[i].im);
        }
    }
}

inline void ifft_odd(uint64_t n, uint64_t p, qd_complex *x, qd_complex *y, qd w[]) {
    uint64_t n2 = n >> 1;
    uint64_t n4 = n >> 2;
    uint64_t l  = n2;
    uint64_t m  = 1;
    uint64_t lp = p - 1;
    uint64_t mp = 0;
#pragma omp target data map(tofrom : x[ : n]) map(to : w[ : n / 4 + 1]) map(alloc : y[ : n])
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
                    inv_butterfly(x[k + (j << mp)],
                                  x[k + (j << mp) + n2],
                                  y[k + (j << (mp + 1))],
                                  y[k + (j << (mp + 1)) + m],
                                  a, b);
                }
            }
            swap(&x, &y);
            l >>= 1;
            m <<= 1;
            mp++;
            lp--;
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x[i].re, n, y[i].re);
            div_pwr2(x[i].im, n, y[i].im);
        }
    }
}
}  // namespace StockhamGPUN4Complex

void stockham_gpu_n4_complex(uint64_t n, uint64_t p, qd_complex *x[], qd w[]) {
    qd_complex *y = (qd_complex *)calloc(n, sizeof(qd_complex));
    if (p & 1)
        StockhamGPUN4Complex::fft_odd(n, p, *x, y, w);
    else
        StockhamGPUN4Complex::fft_even(n, p, *x, y, w);
    free(y);
}

void stockham_gpu_n4_complex(uint64_t n, uint64_t p, qd_complex *x[], qd w[], Timer &timer, Timer &h2d_timer, Timer &d2h_timer, Timer &kernel_timer) {
    qd_complex *y = (qd_complex *)calloc(n, sizeof(qd_complex));
    timer.start();
    if (p & 1)
        StockhamGPUN4Complex::fft_odd(n, p, *x, y, w, h2d_timer, d2h_timer, kernel_timer);
    else
        StockhamGPUN4Complex::fft_even(n, p, *x, y, w, h2d_timer, d2h_timer, kernel_timer);
    timer.stop();
    free(y);
}

void inv_stockham_gpu_n4_complex(uint64_t n, uint64_t p, qd_complex *x[], qd w[]) {
    qd_complex *y = (qd_complex *)calloc(n, sizeof(qd_complex));
    if (p & 1)
        StockhamGPUN4Complex::ifft_odd(n, p, *x, y, w);
    else
        StockhamGPUN4Complex::ifft_even(n, p, *x, y, w);
    free(y);
}