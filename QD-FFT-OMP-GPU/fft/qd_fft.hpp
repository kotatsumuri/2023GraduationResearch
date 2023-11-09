#pragma once
#include <cstdint>

#include "../qd/qd.hpp"
#include "fft_util.hpp"

void stockham(uint32_t n, uint32_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint32_t l = n >> 1;
    uint32_t m = 1;

    for (uint32_t t = 0; t < p; t++) {
        #pragma omp target teams distribute parallel for
        for (uint32_t j = 0; j < l; j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for (uint32_t k = 0; k < m; k++) {
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

void inv_stockham(uint32_t n, uint32_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint32_t l = n >> 1;
    uint32_t m = 1;

    for (uint32_t t = 0; t < p; t++) {
        for (uint32_t j = 0; j < l; j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for (uint32_t k = 0; k < m; k++) {
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

    for (uint32_t i = 0; i < n; i++) {
        div_pwr2(x0[i], n, x0[i]);
        div_pwr2(ix0[i], n, ix0[i]);
    }
}

void dif(uint32_t N, uint32_t n, double *x[], double *ix[], double *y[], double *iy[], qd cos_table[], qd sin_table[]) {
    if (n <= 1)
        return;

    for (uint32_t i = 0; i < (n >> 1); i++) {
        uint32_t j = i + (n >> 1);
        add(x[i], x[j], y[i]);
        add(ix[i], ix[j], iy[i]);

        sub(x[i], x[j], y[j]);
        sub(ix[i], ix[j], iy[j]);

        double *a = (double *)cos_table[N / n * i];
        double *b = (double *)sin_table[N / n * i];

        mul(y[j], a, x[i]);
        mul(y[j], b, ix[i]);
        mul(iy[j], a, x[j]);
        mul(iy[j], b, ix[j]);
        add(x[i], ix[j], y[j]);
        sub(x[j], ix[i], iy[j]);
    }

    dif(N, n >> 1, y, iy, x, ix, cos_table, sin_table);
    dif(N, n >> 1, y + (n >> 1), iy + (n >> 1), x + (n >> 1), ix + (n >> 1), cos_table, sin_table);

    for (uint32_t i = 0; i < (n >> 1); i++) {
        swap(y + i, x + 2 * i);
        swap(iy + i, ix + 2 * i);
        swap(y + i + (n >> 1), x + 2 * i + 1);
        swap(iy + i + (n >> 1), ix + 2 * i + 1);
    }
}

void inv_dif(uint32_t N, uint32_t n, double *x[], double *ix[], double *y[], double *iy[], qd cos_table[], qd sin_table[]) {
    if (n <= 1) {
        div_pwr2(x[0], N, y[0]);
        div_pwr2(ix[0], N, iy[0]);
        return;
    }

    for (uint32_t i = 0; i < (n >> 1); i++) {
        uint32_t j = i + (n >> 1);
        add(x[i], x[j], y[i]);
        add(ix[i], ix[j], iy[i]);

        sub(x[i], x[j], y[j]);
        sub(ix[i], ix[j], iy[j]);

        double *a = (double *)cos_table[N / n * i];
        double *b = (double *)sin_table[N / n * i];

        mul(y[j], a, x[i]);
        mul(y[j], b, ix[i]);
        mul(iy[j], a, x[j]);
        mul(iy[j], b, ix[j]);
        sub(x[i], ix[j], y[j]);
        add(x[j], ix[i], iy[j]);
    }

    inv_dif(N, n >> 1, y, iy, x, ix, cos_table, sin_table);
    inv_dif(N, n >> 1, y + (n >> 1), iy + (n >> 1), x + (n >> 1), ix + (n >> 1),
            cos_table, sin_table);

    for (uint32_t i = 0; i < (n >> 1); i++) {
        swap(y + i, x + 2 * i);
        swap(iy + i, ix + 2 * i);
        swap(y + i + (n >> 1), x + 2 * i + 1);
        swap(iy + i + (n >> 1), ix + 2 * i + 1);
    }
}

inline void stockham_for_sixstep(uint32_t n, uint32_t p, uint32_t u, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint32_t l = n >> 1;
    uint32_t m = 1;

    for (uint32_t t = 0; t < p; t++) {
        for (uint32_t j = 0; j < l; j++) {
            double *a = (double *)cos_table[j * n / (2 * l) * u];
            double *b = (double *)sin_table[j * n / (2 * l) * u];
            for (uint32_t k = 0; k < m; k++) {
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

    if ((p % 2)) {
        for (uint32_t i = 0;i < n;i++) {
            copy(y[i], x[i]);
            copy(iy[i], ix[i]);
        }
    }
}

void sixstep_fft(uint32_t N, uint32_t logN, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]) {
    uint32_t n1 = 1 << (logN / 2);
    uint32_t n2 = 1 << ((logN + 1) / 2);
    uint32_t u1 = N / n1;
    uint32_t u2 = N / n2;
    uint32_t harf_logN1 = logN / 2;
    uint32_t harf_logN2 = (logN + 1) / 2;

    #pragma omp target data map(to:x[:N], ix[:N], cos_table[:N], sin_table[:N]) map(from:y[:N], iy[:N])
    {
        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t i = 0;i < 1;i++) {
            for (uint32_t j = 0;j < 1;j++) {
                copy(x[i * n2 + j], y[j * n1 + i]);
                copy(ix[i * n2 + j], iy[j * n1 + i]);
            }
        }

        #pragma omp target teams distribute parallel for
        for (uint32_t j = 0;j < n2;j++) {
            stockham_for_sixstep(n1, harf_logN1, u1, y + j * n1, iy + j * n1, x + j * n1, ix + j * n1, cos_table, sin_table);
        }

        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t j = 0;j < n2;j++) {
            for (uint32_t i = 0;i < n1;i++) {
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
        for (uint32_t i = 0;i < n1;i++) {
            stockham_for_sixstep(n2, harf_logN2, u2, x + i * n2, ix + i * n2, y + i * n2, iy + i * n2, cos_table, sin_table);
        }

        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t i = 0;i < n1;i++) {
            for (uint32_t j = 0;j < n2;j++) {
                copy(x[i * n2 + j], y[j * n1 + i]);
                copy(ix[i * n2 + j], iy[j * n1 + i]);
            }
        }
    }
}

inline void stockham(uint32_t n, uint32_t p, double x[], double y[], double cos_table[], double sin_table[]) {
    double *x0     = x;
    double *x1     = y;
    uint32_t l = n >> 1;
    uint32_t m = 1;

    for (uint32_t t = 0; t < p; t++) {
        for (uint32_t j = 0; j < l; j++) {
            double *a = cos_table + j * n / (2 * l) * 4;
            double *b = sin_table + j * n / (2 * l) * 4;
            for (uint32_t k = 0; k < m; k++) {
                double *c0  = x0 + (k + j * m) * 8 ;
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
        for (uint32_t i = 0;i < n;i++) {
            copy(y + i * 8, x + i * 8);
            copy(y + i * 8 + 4, x + i * 8 + 4);
        }
    }
}

inline void inv_stockham(uint32_t n, uint32_t p, double x[], double y[], double cos_table[], double sin_table[]) {
    double *x0     = x;
    double *x1     = y;
    uint32_t l = n >> 1;
    uint32_t m = 1;

    for (uint32_t t = 0; t < p; t++) {
        for (uint32_t j = 0; j < l; j++) {
            double *a = cos_table + j * n / (2 * l) * 4;
            double *b = sin_table + j * n / (2 * l) * 4;
            for (uint32_t k = 0; k < m; k++) {
                double *c0  = x0 + (k + j * m) * 8 ;
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
        for (uint32_t i = 0;i < n;i++) {
            div_pwr2(y + i * 8, n, x + i * 8);
            div_pwr2(y + i * 8 + 4, n, x + i * 8 + 4);
        }
    }
    else {
        for (uint32_t i = 0;i < n;i++) {
            div_pwr2(x + i * 8, n, x + i * 8);
            div_pwr2(x + i * 8 + 4, n, x + i * 8 + 4);
        }
    }
}

void sixstep_fft(uint32_t N, uint32_t logN, double x[],  double y[], double cos_table_N[], double sin_table_N[], double cos_table_n1[], double sin_table_n1[], double cos_table_n2[], double sin_table_n2[]) {
    uint32_t n1 = 1 << (logN / 2);
    uint32_t n2 = 1 << ((logN + 1) / 2);
    uint32_t harf_logN1 = logN / 2;
    uint32_t harf_logN2 = (logN + 1) / 2;

    #pragma omp target data map(to:x[:8 * N], cos_table_N[:4 * N], sin_table_N[:4 * N], cos_table_n1[:4 * n1], sin_table_n1[:4 * n1], cos_table_n2[:4 * n2], sin_table_n2[:4 * n2]) map(from:y[:8 * N])
    {
        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t i = 0;i < n1;i++) {
            for (uint32_t j = 0;j < n2;j++) {
                copy(x + (i * n2 + j) * 8, y + (j * n1 + i) * 8);
                copy(x + (i * n2 + j) * 8 + 4, y + (j * n1 + i) * 8 + 4);
            }
        }

        #pragma omp target teams distribute parallel for
        for (uint32_t j = 0;j < n2;j++) {
            stockham(n1, harf_logN1, y + j * n1 * 8, x + j * n1 * 8, cos_table_n1, sin_table_n1);
        }

        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t j = 0;j < n2;j++) {
            for (uint32_t i = 0;i < n1;i++) {
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
        for (uint32_t i = 0;i < n1;i++) {
            stockham(n2, harf_logN2, x + i * n2 * 8, y + i * n2 * 8, cos_table_n2, sin_table_n2);
        }

        #pragma omp target teams distribute parallel for collapse(2)
        for (uint32_t i = 0;i < n1;i++) {
            for (uint32_t j = 0;j < n2;j++) {
                copy(x + (i * n2 + j) * 8, y + (j * n1 + i) * 8);
                copy(x + (i * n2 + j) * 8 + 4, y + (j * n1 + i) * 8 + 4);
            }
        }
    }
}

