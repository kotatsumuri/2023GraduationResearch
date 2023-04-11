#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include "Utility.h"

void testMakeW() {
    int n = 16;
    double complex* w = makeW(n);
    int i;
    for(i = 0;i < 16;i++) {
        assert(creal(w[i]) == cos(2.0 * M_PI * (double)i / (double)n));
        assert(cimag(w[i]) == - sin(2.0 * M_PI * (double)i / (double)n));
        printComplex(w[i]);
        puts("");
    }
}

void testBitReversed() {
    int n = 16;
    int i;
    double complex* x = (double complex *)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++) {
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)(n * 2));
        printComplex(x[i]);
        puts("");
    }
    bitReversed(n, x);
    for(i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }
}

int main() {
    testMakeW();
    testBitReversed();
    return 0;
}