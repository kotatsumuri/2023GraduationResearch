#pragma once
#include <cstdint>

#include "../qd/qd.hpp"
#include "fft_util.hpp"

#define QD_FFT_IN_SWAP

void stockham(uint64_t n, uint64_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
#pragma omp target teams distribute parallel for
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for (uint64_t k = 0; k < m; k++) {
                double *c0  = (double *)x0[k + j * m];
                double *ic0 = (double *)ix0[k + j * m];

                double *c1  = (double *)x0[k + j * m + l * m];
                double *ic1 = (double *)ix0[k + j * m + l * m];

                double *xt1  = (double *)x1[k + 2 * j * m + m];
                double *ixt1 = (double *)ix1[k + 2 * j * m + m];

                // c0 + c1
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1[k + 2 * j * m]);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, ix1[k + 2 * j * m]);

                // (a - bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) + b * (ic0 - ic1) + j * (a * (ic0 - ic1) - b
                // * (c0 - c1)) xt1 <- c0 - c1
                sub(c0, c1, xt1);
                // ixt1 <- ic0 - ic1
                sub(ic0, ic1, ixt1);
                // c0 <- a * xt1 = a * (c0 - c1)
                mul(a, xt1, c0);
                // ic0 <- a * ixt1 = a * (ic0 - ic1)
                mul(a, ixt1, ic0);
                // c1 <- b * ixt1 = b * (ic0 - ic1)
                mul(b, ixt1, c1);
                // ic1 <- b * xt1 = b * (c0 - c1)
                mul(b, xt1, ic1);
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) + b * (ic0 -
                // ic1)
                add(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) - b *
                // (c0 - c1)
                sub(ic0, ic1, ixt1);
            }
        }
        swap(&x0, &x1);
        swap(&ix0, &ix1);

        l >>= 1;
        m <<= 1;
    }
}

void inv_stockham(uint64_t n, uint64_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for (uint64_t k = 0; k < m; k++) {
                double *c0  = (double *)x0[k + j * m];
                double *ic0 = (double *)ix0[k + j * m];

                double *c1  = (double *)x0[k + j * m + l * m];
                double *ic1 = (double *)ix0[k + j * m + l * m];

                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];

                // c0 + c1
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1[k + 2 * j * m]);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, ix1[k + 2 * j * m]);

                // (a + bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) - b * (ic0 - ic1) + j * (a * (ic0 - ic1) + b
                // * (c0 - c1)) xt1 <- c0 - c1
                sub(c0, c1, xt1);
                // ixt1 <- ic0 - ic1
                sub(ic0, ic1, ixt1);
                // c0 <- a * xt1 = a * (c0 - c1)
                mul(a, xt1, c0);
                // ic0 <- a * ixt1 = a * (ic0 - ic1)
                mul(a, ixt1, ic0);
                // c1 <- b * ixt1 = b * (ic0 - ic1)
                mul(b, ixt1, c1);
                // ic1 <- b * xt1 = b * (c0 - c1)
                mul(b, xt1, ic1);
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) - b * (ic0 -
                // ic1)
                sub(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) + b *
                // (c0 - c1)
                add(ic0, ic1, ixt1);
            }
        }
        swap(&x0, &x1);
        swap(&ix0, &ix1);

        l >>= 1;
        m <<= 1;
    }

    for (uint64_t i = 0; i < n; i++) {
        div_pwr2(x0[i], n, x0[i]);
        div_pwr2(ix0[i], n, ix0[i]);
    }
}

inline void stockham_for_sixstep(uint64_t n, uint64_t p, uint64_t u, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0           = x;
    qd *ix0          = ix;
    qd *x1           = y;
    qd *ix1          = iy;
    uint64_t n_harf  = n >> 1ull;
    uint64_t l       = n_harf;
    uint64_t m       = 1ull;
    uint64_t ab_diff = u;

    for (uint64_t t = 0; t < p; t++) {
        uint64_t ab_i  = 0;
        uint64_t c0_i  = 0;
        uint64_t c1_i  = n_harf;
        uint64_t xt1_i = m;
        uint64_t x1_i  = 0;
        for (uint64_t j = 0; j < l; j++) {
            // double *a = (double *)cos_table[j * n / (2 * l) * u];
            // double *b = (double *)sin_table[j * n / (2 * l) * u];
            double *a = (double *)cos_table[ab_i];
            double *b = (double *)sin_table[ab_i];
            for (uint64_t k = 0; k < m; k++) {
                // double *c0  = (double *)x0[k + j * m];
                // double *ic0 = (double *)ix0[k + j * m];
                double *c0  = (double *)x0[c0_i];
                double *ic0 = (double *)ix0[c0_i];

                // double *c1  = (double *)x0[k + j * m + l * m];
                // double *ic1 = (double *)ix0[k + j * m + l * m];
                double *c1  = (double *)x0[c1_i];
                double *ic1 = (double *)ix0[c1_i];

                // double *xt1  = (double *)x1[k + 2 * j * m + m];
                // double *ixt1 = (double *)ix1[k + 2 * j * m + m];
                double *xt1  = (double *)x1[xt1_i];
                double *ixt1 = (double *)ix1[xt1_i];

                // // c0 + c1
                // // x1[k + 2jm] <- c0 + c1
                // add(c0, c1, x1[k + 2 * j * m]);
                // // ix1[k + 2jm] <- ic0 + ic1
                // add(ic0, ic1, ix1[k + 2 * j * m]);
                // c0 + c1
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1[x1_i]);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, ix1[x1_i]);

                // (a - bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) + b * (ic0 - ic1) + j * (a * (ic0 - ic1) - b
                // * (c0 - c1)) xt1 <- c0 - c1
                sub(c0, c1, xt1);
                // ixt1 <- ic0 - ic1
                sub(ic0, ic1, ixt1);
                // c0 <- a * xt1 = a * (c0 - c1)
                mul(a, xt1, c0);
                // ic0 <- a * ixt1 = a * (ic0 - ic1)
                mul(a, ixt1, ic0);
                // c1 <- b * ixt1 = b * (ic0 - ic1)
                mul(b, ixt1, c1);
                // ic1 <- b * xt1 = b * (c0 - c1)
                mul(b, xt1, ic1);
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) + b * (ic0 -
                // ic1)
                add(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) - b *
                // (c0 - c1)
                sub(ic0, ic1, ixt1);

                c0_i++;
                c1_i++;
                xt1_i++;
                x1_i++;
            }
            ab_i += ab_diff;
            xt1_i += m;
            x1_i += m;
        }
        swap(&x0, &x1);
        swap(&ix0, &ix1);

        l >>= 1ull;
        m <<= 1ull;
        ab_diff <<= 1ull;
    }
#ifdef QD_FFT_IN_SWAP
    if ((p % 2)) {
        for (uint64_t i = 0; i < n; i++) {
            copy(y[i], x[i]);
            copy(iy[i], ix[i]);
        }
    }
#endif
}

void sixstep_fft(uint64_t N, uint64_t logN, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    uint64_t harf_logN1 = logN >> 1ull;
    uint64_t harf_logN2 = (logN + 1ull) >> 1ull;
    uint64_t n1         = 1ull << harf_logN1;
    uint64_t n2         = 1ull << harf_logN2;
    uint64_t u1         = N >> harf_logN1;
    uint64_t u2         = N >> harf_logN2;

#pragma omp target data map(to : x[ : N], ix[ : N], cos_table[ : N], sin_table[ : N]) map(from : y[ : N], iy[ : N])
    {
#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t i = 0; i < n1; i++) {
            for (uint64_t j = 0; j < n2; j++) {
                copy(x[i * n2 + j], y[j * n1 + i]);
                copy(ix[i * n2 + j], iy[j * n1 + i]);
            }
        }

#pragma omp target teams distribute parallel for
        for (uint64_t j = 0; j < n2; j++) {
            stockham_for_sixstep(n1, harf_logN1, u1, y + j * n1, iy + j * n1, x + j * n1, ix + j * n1, cos_table, sin_table);
        }
#ifndef QD_FFT_NO_SWAP
#ifdef QD_FFT_OUT_SWAP
        if (harf_logN1 & 1) {
#pragma omp target teams distribute parallel for
            for (uint64_t i = 0; i < N; i++) {
                copy(x[i], y[i]);
                copy(ix[i], iy[i]);
            }
        }
#endif

#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t j = 0; j < n2; j++) {
            for (uint64_t i = 0; i < n1; i++) {
                qd tmp0, tmp1;
                double *a = (double *)cos_table[i * j];
                double *b = (double *)sin_table[i * j];
                mul(y[j * n1 + i], a, x[i * n2 + j]);
                mul(y[j * n1 + i], b, tmp0);
                mul(iy[j * n1 + i], a, ix[i * n2 + j]);
                mul(iy[j * n1 + i], b, tmp1);
                add(x[i * n2 + j], tmp1, x[i * n2 + j]);
                sub(ix[i * n2 + j], tmp0, ix[i * n2 + j]);
            }
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n1; i++) {
            stockham_for_sixstep(n2, harf_logN2, u2, x + i * n2, ix + i * n2, y + i * n2, iy + i * n2, cos_table, sin_table);
        }

#ifdef QD_FFT_OUT_SWAP
        if (harf_logN2 & 1) {
#pragma omp target teams distribute parallel for
            for (uint64_t i = 0; i < N; i++) {
                copy(y[i], x[i]);
                copy(iy[i], ix[i]);
            }
        }
#endif

#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t i = 0; i < n1; i++) {
            for (uint64_t j = 0; j < n2; j++) {
                copy(x[i * n2 + j], y[j * n1 + i]);
                copy(ix[i * n2 + j], iy[j * n1 + i]);
            }
        }
#else
        if (harf_logN1 & 1) {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < n2; j++) {
                for (uint64_t i = 0; i < n1; i++) {
                    qd tmp0, tmp1;
                    double *a = (double *)cos_table[i * j];
                    double *b = (double *)sin_table[i * j];
                    mul(x[j * n1 + i], a, y[i * n2 + j]);
                    mul(x[j * n1 + i], b, tmp0);
                    mul(ix[j * n1 + i], a, iy[i * n2 + j]);
                    mul(ix[j * n1 + i], b, tmp1);
                    add(y[i * n2 + j], tmp1, y[i * n2 + j]);
                    sub(iy[i * n2 + j], tmp0, iy[i * n2 + j]);
                }
            }

#pragma omp target teams distribute parallel for
            for (uint64_t i = 0; i < n1; i++) {
                stockham_for_sixstep(n2, harf_logN2, u2, y + i * n2, iy + i * n2, x + i * n2, ix + i * n2, cos_table, sin_table);
            }

            if (harf_logN2 & 1) {
#pragma omp target teams distribute parallel for collapse(2)
                for (uint64_t i = 0; i < n1; i++) {
                    for (uint64_t j = 0; j < n2; j++) {
                        copy(x[i * n2 + j], y[j * n1 + i]);
                        copy(ix[i * n2 + j], iy[j * n1 + i]);
                    }
                }
            } else {
#pragma omp target teams distribute parallel for collapse(2)
                for (uint64_t i = 0; i < n1; i++) {
                    for (uint64_t j = 0; j < n2; j++) {
                        copy(y[i * n2 + j], x[j * n1 + i]);
                        copy(iy[i * n2 + j], ix[j * n1 + i]);
                    }
                }
            }
        } else {
#pragma omp target teams distribute parallel for collapse(2)
            for (uint64_t j = 0; j < n2; j++) {
                for (uint64_t i = 0; i < n1; i++) {
                    qd tmp0, tmp1;
                    double *a = (double *)cos_table[i * j];
                    double *b = (double *)sin_table[i * j];
                    mul(y[j * n1 + i], a, x[i * n2 + j]);
                    mul(y[j * n1 + i], b, tmp0);
                    mul(iy[j * n1 + i], a, ix[i * n2 + j]);
                    mul(iy[j * n1 + i], b, tmp1);
                    add(x[i * n2 + j], tmp1, x[i * n2 + j]);
                    sub(ix[i * n2 + j], tmp0, ix[i * n2 + j]);
                }
            }

#pragma omp target teams distribute parallel for
            for (uint64_t i = 0; i < n1; i++) {
                stockham_for_sixstep(n2, harf_logN2, u2, x + i * n2, ix + i * n2, y + i * n2, iy + i * n2, cos_table, sin_table);
            }

            if (harf_logN2 & 1) {
#pragma omp target teams distribute parallel for collapse(2)
                for (uint64_t i = 0; i < n1; i++) {
                    for (uint64_t j = 0; j < n2; j++) {
                        copy(y[i * n2 + j], x[j * n1 + i]);
                        copy(iy[i * n2 + j], ix[j * n1 + i]);
                    }
                }
            } else {
#pragma omp target teams distribute parallel for collapse(2)
                for (uint64_t i = 0; i < n1; i++) {
                    for (uint64_t j = 0; j < n2; j++) {
                        copy(x[i * n2 + j], y[j * n1 + i]);
                        copy(ix[i * n2 + j], iy[j * n1 + i]);
                    }
                }
            }
        }
#endif
    }
}

inline void stockham(uint64_t n, uint64_t p, double x[], double y[], double cos_table[], double sin_table[]) {
    double *x0 = x;
    double *x1 = y;
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = cos_table + j * n / (2 * l) * 4;
            double *b = sin_table + j * n / (2 * l) * 4;
            for (uint64_t k = 0; k < m; k++) {
                double *c0  = x0 + (k + j * m) * 8;
                double *ic0 = (double *)x0 + (k + j * m) * 8 + 4;

                double *c1  = (double *)x0 + (k + j * m + l * m) * 8;
                double *ic1 = (double *)x0 + (k + j * m + l * m) * 8 + 4;

                double *xt1  = (double *)x1 + (k + 2 * j * m + m) * 8;
                double *ixt1 = (double *)x1 + (k + 2 * j * m + m) * 8 + 4;

                // c0 + c1
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1 + (k + 2 * j * m) * 8);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, x1 + (k + 2 * j * m) * 8 + 4);

                // (a - bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) + b * (ic0 - ic1) + j * (a * (ic0 - ic1) - b
                // * (c0 - c1)) xt1 <- c0 - c1
                sub(c0, c1, xt1);
                // ixt1 <- ic0 - ic1
                sub(ic0, ic1, ixt1);
                // c0 <- a * xt1 = a * (c0 - c1)
                mul(a, xt1, c0);
                // ic0 <- a * ixt1 = a * (ic0 - ic1)
                mul(a, ixt1, ic0);
                // c1 <- b * ixt1 = b * (ic0 - ic1)
                mul(b, ixt1, c1);
                // ic1 <- b * xt1 = b * (c0 - c1)
                mul(b, xt1, ic1);
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) + b * (ic0 -
                // ic1)
                add(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) - b *
                // (c0 - c1)
                sub(ic0, ic1, ixt1);
            }
        }
        swap(&x0, &x1);

        l >>= 1;
        m <<= 1;
    }

    if ((p % 2)) {
        for (uint64_t i = 0; i < n; i++) {
            copy(y + i * 8, x + i * 8);
            copy(y + i * 8 + 4, x + i * 8 + 4);
        }
    }
}

inline void inv_stockham(uint64_t n, uint64_t p, double x[], double y[], double cos_table[], double sin_table[]) {
    double *x0 = x;
    double *x1 = y;
    uint64_t l = n >> 1;
    uint64_t m = 1;

    for (uint64_t t = 0; t < p; t++) {
        for (uint64_t j = 0; j < l; j++) {
            double *a = cos_table + j * n / (2 * l) * 4;
            double *b = sin_table + j * n / (2 * l) * 4;
            for (uint64_t k = 0; k < m; k++) {
                double *c0  = x0 + (k + j * m) * 8;
                double *ic0 = (double *)x0 + (k + j * m) * 8 + 4;

                double *c1  = (double *)x0 + (k + j * m + l * m) * 8;
                double *ic1 = (double *)x0 + (k + j * m + l * m) * 8 + 4;

                double *xt1  = (double *)x1 + (k + 2 * j * m + m) * 8;
                double *ixt1 = (double *)x1 + (k + 2 * j * m + m) * 8 + 4;

                // c0 + c1
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1 + (k + 2 * j * m) * 8);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, x1 + (k + 2 * j * m) * 8 + 4);

                // (a + bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) - b * (ic0 - ic1) + j * (a * (ic0 - ic1) + b
                // * (c0 - c1)) xt1 <- c0 - c1
                sub(c0, c1, xt1);
                // ixt1 <- ic0 - ic1
                sub(ic0, ic1, ixt1);
                // c0 <- a * xt1 = a * (c0 - c1)
                mul(a, xt1, c0);
                // ic0 <- a * ixt1 = a * (ic0 - ic1)
                mul(a, ixt1, ic0);
                // c1 <- b * ixt1 = b * (ic0 - ic1)
                mul(b, ixt1, c1);
                // ic1 <- b * xt1 = b * (c0 - c1)
                mul(b, xt1, ic1);
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) - b * (ic0 -
                // ic1)
                sub(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) + b *
                // (c0 - c1)
                add(ic0, ic1, ixt1);
            }
        }
        swap(&x0, &x1);

        l >>= 1;
        m <<= 1;
    }

    if ((p % 2)) {
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(y + i * 8, n, x + i * 8);
            div_pwr2(y + i * 8 + 4, n, x + i * 8 + 4);
        }
    } else {
        for (uint64_t i = 0; i < n; i++) {
            div_pwr2(x + i * 8, n, x + i * 8);
            div_pwr2(x + i * 8 + 4, n, x + i * 8 + 4);
        }
    }
}

void sixstep_fft(uint64_t N, uint64_t logN, double x[], double y[], double cos_table_N[], double sin_table_N[], double cos_table_n1[], double sin_table_n1[], double cos_table_n2[], double sin_table_n2[]) {
    uint64_t n1         = 1 << (logN / 2);
    uint64_t n2         = 1 << ((logN + 1) / 2);
    uint64_t harf_logN1 = logN / 2;
    uint64_t harf_logN2 = (logN + 1) / 2;

#pragma omp target data map(to : x[ : 8 * N], cos_table_N[ : 4 * N], sin_table_N[ : 4 * N], cos_table_n1[ : 4 * n1], sin_table_n1[ : 4 * n1], cos_table_n2[ : 4 * n2], sin_table_n2[ : 4 * n2]) map(from : y[ : 8 * N])
    {
#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t i = 0; i < n1; i++) {
            for (uint64_t j = 0; j < n2; j++) {
                copy(x + (i * n2 + j) * 8, y + (j * n1 + i) * 8);
                copy(x + (i * n2 + j) * 8 + 4, y + (j * n1 + i) * 8 + 4);
            }
        }

#pragma omp target teams distribute parallel for
        for (uint64_t j = 0; j < n2; j++) {
            stockham(n1, harf_logN1, y + j * n1 * 8, x + j * n1 * 8, cos_table_n1, sin_table_n1);
        }

#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t j = 0; j < n2; j++) {
            for (uint64_t i = 0; i < n1; i++) {
                qd tmp0, tmp1;
                double *a = cos_table_N + i * j * 4;
                double *b = sin_table_N + i * j * 4;
                mul(y + (j * n1 + i) * 8, a, x + (i * n2 + j) * 8);
                mul(y + (j * n1 + i) * 8, b, tmp0);
                mul(y + (j * n1 + i) * 8 + 4, a, x + (i * n2 + j) * 8 + 4);
                mul(y + (j * n1 + i) * 8 + 4, b, tmp1);
                add(x + (i * n2 + j) * 8, tmp1, x + (i * n2 + j) * 8);
                sub(x + (i * n2 + j) * 8 + 4, tmp0, x + (i * n2 + j) * 8 + 4);
            }
        }

#pragma omp target teams distribute parallel for
        for (uint64_t i = 0; i < n1; i++) {
            stockham(n2, harf_logN2, x + i * n2 * 8, y + i * n2 * 8, cos_table_n2, sin_table_n2);
        }

#pragma omp target teams distribute parallel for collapse(2)
        for (uint64_t i = 0; i < n1; i++) {
            for (uint64_t j = 0; j < n2; j++) {
                copy(x + (i * n2 + j) * 8, y + (j * n1 + i) * 8);
                copy(x + (i * n2 + j) * 8 + 4, y + (j * n1 + i) * 8 + 4);
            }
        }
    }
}