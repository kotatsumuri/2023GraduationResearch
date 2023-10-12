#pragma once
#include <cstdint>

#include "../qd/qd.hpp"

void stockham(uint16_t n, uint16_t p, qd x[], qd ix[], qd y[], qd iy[],
              qd cos_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint16_t l = n >> 1;
    uint16_t m = 1;

    for (uint16_t t = 0; t < p - 2; t++) {
        uint16_t j       = 0;
        uint16_t cos_ptr = 0;
        uint16_t sin_ptr = n >> 2;
        for (; j < l / 4; j++) {
            double *a = (double *)cos_table[cos_ptr];
            double *b = (double *)cos_table[sin_ptr];
            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                add(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
            }
            cos_ptr += m;
            sin_ptr -= m;
        }

        for (; j < l / 2; j++) {
            qd a;
            minus(cos_table[cos_ptr], a);
            double *b = (double *)cos_table[sin_ptr];

            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                add(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
            }
            cos_ptr -= m;
            sin_ptr += m;
        }

        for (; j < 3 * l / 4; j++) {
            qd a;
            minus(cos_table[cos_ptr], a);
            qd b;
            minus(cos_table[sin_ptr], b);
            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                add(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
            }
            cos_ptr += m;
            sin_ptr -= m;
        }

        for (; j < l; j++) {
            double *a = (double *)cos_table[cos_ptr];
            qd b;
            minus(cos_table[sin_ptr], b);

            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                add(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
            }
            cos_ptr -= m;
            sin_ptr += m;
        }
        qd *tmp = x0;
        x0      = x1;
        x1      = tmp;
        tmp     = ix0;
        ix0     = ix1;
        ix1     = tmp;

        l >>= 1;
        m <<= 1;
    }

    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k];
        double *ic0  = (double *)ix0[k];
        double *c1   = (double *)x0[k + (m << 1)];
        double *ic1  = (double *)ix0[k + (m << 1)];
        double *xt1  = (double *)x1[k + m];
        double *ixt1 = (double *)ix1[k + m];
        add(c0, c1, x1[k]);
        add(ic0, ic1, ix1[k]);
        sub(c0, c1, xt1);
        sub(ic0, ic1, ixt1);
    }

    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k + m];
        double *ic0  = (double *)ix0[k + m];
        double *c1   = (double *)x0[k + 3 * m];
        double *ic1  = (double *)ix0[k + 3 * m];
        double *xt1  = (double *)x1[k + 3 * m];
        double *ixt1 = (double *)ix1[k + 3 * m];
        add(c0, c1, x1[k + (m << 1)]);
        add(ic0, ic1, ix1[k + (m << 1)]);
        sub(c1, c0, ixt1);
        sub(ic0, ic1, xt1);
    }
    qd *tmp = x0;
    x0      = x1;
    x1      = tmp;
    tmp     = ix0;
    ix0     = ix1;
    ix1     = tmp;

    m <<= 1;
    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k];
        double *ic0  = (double *)ix0[k];
        double *c1   = (double *)x0[k + m];
        double *ic1  = (double *)ix0[k + m];
        double *xt1  = (double *)x1[k + m];
        double *ixt1 = (double *)ix1[k + m];
        add(c0, c1, x1[k]);
        add(ic0, ic1, ix1[k]);
        sub(c0, c1, xt1);
        sub(ic0, ic1, ixt1);
    }
    tmp = x0;
    x0  = x1;
    x1  = tmp;
    tmp = ix0;
    ix0 = ix1;
    ix1 = tmp;
}

void inv_stockham(uint16_t n, uint16_t p, qd x[], qd ix[], qd y[], qd iy[],
                  qd cos_table[]) {
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint16_t l = n >> 1;
    uint16_t m = 1;

    for (uint16_t t = 0; t < p - 2; t++) {
        uint16_t j       = 0;
        uint16_t cos_ptr = 0;
        uint16_t sin_ptr = n >> 2;
        for (; j < l / 4; j++) {
            double *a = (double *)cos_table[cos_ptr];
            double *b = (double *)cos_table[sin_ptr];
            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                sub(c0, c1, xt1);
                add(ic0, ic1, ixt1);
            }
            cos_ptr += m;
            sin_ptr -= m;
        }

        for (; j < l / 2; j++) {
            qd a;
            minus(cos_table[cos_ptr], a);
            double *b = (double *)cos_table[sin_ptr];

            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                sub(c0, c1, xt1);
                add(ic0, ic1, ixt1);
            }
            cos_ptr -= m;
            sin_ptr += m;
        }

        for (; j < 3 * l / 4; j++) {
            qd a;
            minus(cos_table[cos_ptr], a);
            qd b;
            minus(cos_table[sin_ptr], b);
            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                sub(c0, c1, xt1);
                add(ic0, ic1, ixt1);
            }
            cos_ptr += m;
            sin_ptr -= m;
        }

        for (; j < l; j++) {
            double *a = (double *)cos_table[cos_ptr];
            qd b;
            minus(cos_table[sin_ptr], b);

            for (uint16_t k = 0; k < m; k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];
                add(c0, c1, x1[k + ((j * m) << 1)]);
                add(ic0, ic1, ix1[k + ((j * m) << 1)]);
                sub(c0, c1, xt1);
                sub(ic0, ic1, ixt1);
                mul(a, xt1, c0);
                mul(a, ixt1, ic0);
                mul(b, ixt1, c1);
                mul(b, xt1, ic1);
                sub(c0, c1, xt1);
                add(ic0, ic1, ixt1);
            }
            cos_ptr -= m;
            sin_ptr += m;
        }
        qd *tmp = x0;
        x0      = x1;
        x1      = tmp;
        tmp     = ix0;
        ix0     = ix1;
        ix1     = tmp;

        l >>= 1;
        m <<= 1;
    }

    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k];
        double *ic0  = (double *)ix0[k];
        double *c1   = (double *)x0[k + (m << 1)];
        double *ic1  = (double *)ix0[k + (m << 1)];
        double *xt1  = (double *)x1[k + m];
        double *ixt1 = (double *)ix1[k + m];
        add(c0, c1, x1[k]);
        add(ic0, ic1, ix1[k]);
        sub(c0, c1, xt1);
        sub(ic0, ic1, ixt1);
    }

    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k + m];
        double *ic0  = (double *)ix0[k + m];
        double *c1   = (double *)x0[k + 3 * m];
        double *ic1  = (double *)ix0[k + 3 * m];
        double *xt1  = (double *)x1[k + 3 * m];
        double *ixt1 = (double *)ix1[k + 3 * m];
        add(c0, c1, x1[k + (m << 1)]);
        add(ic0, ic1, ix1[k + (m << 1)]);
        sub(c0, c1, ixt1);
        sub(ic1, ic0, xt1);
    }
    qd *tmp = x0;
    x0      = x1;
    x1      = tmp;
    tmp     = ix0;
    ix0     = ix1;
    ix1     = tmp;

    m <<= 1;
    for (int k = 0; k < m; k++) {
        double *c0   = (double *)x0[k];
        double *ic0  = (double *)ix0[k];
        double *c1   = (double *)x0[k + m];
        double *ic1  = (double *)ix0[k + m];
        double *xt1  = (double *)x1[k + m];
        double *ixt1 = (double *)ix1[k + m];
        add(c0, c1, x1[k]);
        add(ic0, ic1, ix1[k]);
        sub(c0, c1, xt1);
        sub(ic0, ic1, ixt1);
    }
    tmp = x0;
    x0  = x1;
    x1  = tmp;
    tmp = ix0;
    ix0 = ix1;
    ix1 = tmp;

    for(int i = 0;i < n;i++) {
        mul_pwr2(x0[i], 1.0 / n, x0[i]);
        mul_pwr2(ix0[i], 1.0 / n, ix0[i]);
    }
}