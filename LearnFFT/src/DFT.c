#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "DFT.h"

void DFT(int n, double complex* x, double complex* y, double complex* w) {
    int i, j;
    for(i = 0;i < n;i++) {
        y[i] = 0;
        for(j = 0;j < n;j++) {
            y[i] += x[j] * w[(i * j) % n];
        }
    }
}

void invDFT(int n, double complex *x, double complex* y, double complex *w) {
    int i, j;
    for(i = 0;i < n;i++) {
        y[i] = 0;
        for(j = 0;j < n;j++) {
            y[i] += x[j] / w[(i * j) % n];
        }
        y[i] /= (double)n;
    }
}