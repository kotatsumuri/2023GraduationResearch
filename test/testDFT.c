#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "Utility.h"
#include "DFT.h"

int main() {
    int i;
    int n = 16;
    double complex* w = makeW(n);
    double complex* x = (double complex*)calloc(n, sizeof(double complex));
    for(i = 0;i < n;i++) {
        x[i] = (double complex)cos(2.0 * M_PI * (double)i / (double)n);
        printComplex(x[i]);
        puts("");
    }

    double complex* y = (double complex*)calloc(n, sizeof(double complex));

    DFT(n, x, y, w);
    
    for(i = 0;i < n;i++) {
        printComplex(y[i]);
        puts("");
    }

    invDFT(n, y, x, w);

    for(i = 0;i < n;i++) {
        printComplex(x[i]);
        puts("");
    }

    return 0;
}