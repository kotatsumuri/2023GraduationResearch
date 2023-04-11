#include "Utility.h"

double complex* makeW(int n) {
    double complex* w = (double complex*)calloc(n, sizeof(double complex));
    for(int i = 0;i < n;i++) 
        w[i] = cexp(I * 2.0 * M_PI_2 * (double)i / (double)n);
    return w;
}