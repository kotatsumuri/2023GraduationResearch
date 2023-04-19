#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "Utility.h"
#include "DFT.h"
#include "FFT.h"

void testDecimationInFrequency() {
    printf("testDecimationInFrequency\n");
    int n = 16;
    int i;
    double complex *w = makeW(n);
    double complex *x = (double complex*)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++)
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)n);
    double complex *y = (double complex*)calloc(n, sizeof(double complex));

    decimationInFrequency(n, n, x, y, w);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }

    invDecimationInFrequency(n, n, x, y, w);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
}

void testCooleyTukey() {
    printf("testCooleyTukey\n");
    int n = 16, m = 0, l = 1, i;
    while(n > l) {
        m++;
        l *= 2;
    }
    double complex *w = makeW(n);
    double complex *x = (double complex*)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++)
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)n);

    cooleyTukey(n, m, x, w);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }

    invCooleyTukey(n, m, x, w);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
}

void testStockham() {
    printf("testStockham\n");
    int n = 16, p = 0, l = 1, i;
    while(n > l) {
        p++;
        l *= 2;
    }
    double complex *w = makeW(n);
    double complex *x = (double complex*)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++)
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)n);
    double complex *y = (double complex*)calloc(n, sizeof(double complex));

    stockham(n, p, x, y, w);
    for(int i = 0;i < n;i++) {
        printComplex(y[i]);
        puts("");
    }

    invStockham(n, p, y, x, w);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
}

void testRealFFT() {
    printf("testRealFFT\n");
    int n = 16, m = 0, l = 1, i;
    while(n > l) {
        m++;
        l *= 2;
    }
    double complex *w = makeW(n);
    double* x = (double*)calloc(n, sizeof(double));
    for(i = 0;i < n;i++)
        x[i] = cos(2.0 * M_PI * (double)i / (double)n);
    double complex *y = (double complex*)calloc(n, sizeof(double complex));
    realFFT(n, m, x, y, w);
    for(int i = 0;i < n;i++) {
        printComplex(y[i]);
        puts("");
    }
    for(i = 0;i < n;i++)
        x[i] = creal(y[i]);
    invRealFFT(n, m, x, y, w);
    for(int i = 0;i < n;i++) {
        printComplex(y[i]);
        puts("");
    }
}

void testSplitRadixFFT() {
    printf("testSplitRadixFFT\n");
    int n = 16, i;
    double complex *w = makeW(n);
    double complex* x = (double complex*)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++)
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)n);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
    splitRadixFFT(n, n, x, w);
    bitReversed(n, x);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
    invSplitRadixFFT(n, n, x, w);
    bitReversed(n, x);
    for(int i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
}

int main() {
    testDecimationInFrequency();
    testCooleyTukey();
    testStockham();
    testRealFFT();
    testSplitRadixFFT();
    return 0;
}