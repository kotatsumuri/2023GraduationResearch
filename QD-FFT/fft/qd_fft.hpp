#pragma once
#include <cstdint>

#include "../qd/qd.hpp"

void stockham(uint16_t n, uint16_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]){
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint16_t l = n >> 1;
    uint16_t m = 1;

    for(uint16_t t = 0;t < p;t++) {
        for(uint16_t j = 0;j < l;j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for(uint16_t k = 0;k < m;k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];
                
                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];
                
                double *xt1  = (double *)x1[k + 2 * j * m + m];
                double *ixt1 = (double *)ix1[k + 2 * j * m + m];
                
                // c0 + c1 
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1[k + 2 * j * m]);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, ix1[k + 2 * j * m]);

                // (a - bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) + b * (ic0 - ic1) + j * (a * (ic0 - ic1) - b * (c0 - c1))
                // xt1 <- c0 - c1
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
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) + b * (ic0 - ic1)
                add(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) - b * (c0 - c1)
                sub(ic0, ic1, ixt1);
            }
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
}

void inv_stockham(uint16_t n, uint16_t p, qd x[], qd ix[], qd y[], qd iy[], qd cos_table[], qd sin_table[]){
    qd *x0     = x;
    qd *ix0    = ix;
    qd *x1     = y;
    qd *ix1    = iy;
    uint16_t l = n >> 1;
    uint16_t m = 1;

    for(uint16_t t = 0;t < p;t++) {
        for(uint16_t j = 0;j < l;j++) {
            double *a = (double *)cos_table[j * n / (2 * l)];
            double *b = (double *)sin_table[j * n / (2 * l)];
            for(uint16_t k = 0;k < m;k++) {
                double *c0   = (double *)x0[k + j * m];
                double *ic0  = (double *)ix0[k + j * m];

                double *c1   = (double *)x0[k + j * m + l * m];
                double *ic1  = (double *)ix0[k + j * m + l * m];

                double *xt1  = (double *)x1[k + ((j * m) << 1) + m];
                double *ixt1 = (double *)ix1[k + ((j * m) << 1) + m];

                // c0 + c1 
                // x1[k + 2jm] <- c0 + c1
                add(c0, c1, x1[k + 2 * j * m]);
                // ix1[k + 2jm] <- ic0 + ic1
                add(ic0, ic1, ix1[k + 2 * j * m]);

                // (a + bj) * (c0 + ic0j - (c1 + ic1j))
                // = a * (c0 - c1) - b * (ic0 - ic1) + j * (a * (ic0 - ic1) + b * (c0 - c1))
                // xt1 <- c0 - c1
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
                // x1[k + 2jm + m] = xt1 <- c0 + c1 = a * (c0 - c1) - b * (ic0 - ic1)
                sub(c0, c1, xt1);
                // ix1[k + 2jm + m] = ixt1 <- ic0 - ic1 = a * (ic0 - ic1) + b * (c0 - c1)
                add(ic0, ic1, ixt1);
            }
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

    for(uint16_t i = 0;i < n;i++) {
        div_pwr2(x0[i], n, x0[i]);
        div_pwr2(ix0[i], n, ix0[i]);
    }
}