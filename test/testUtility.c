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
        assert(creal(w[i]) == cos(2.0 * M_PI_2 * (double)i / (double)n));
        assert(cimag(w[i]) == sin(2.0 * M_PI_2 * (double)i / (double)n));
        printComplex(w[i]);
        puts("");
    }
}

int main() {
    testMakeW();
    return 0;
}