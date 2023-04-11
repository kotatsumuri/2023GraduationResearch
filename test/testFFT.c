#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "Utility.h"
#include "DFT.h"
#include "FFT.h"

void testDecimationInFrequency() {
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

int main() {
    testDecimationInFrequency();
    testCooleyTukey();
    return 0;
}