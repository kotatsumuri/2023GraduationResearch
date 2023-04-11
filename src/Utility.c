#include "Utility.h"

void printComplex(double complex c) {
    if(cimag(c) >= 0)
        printf("%lf + i%lf", creal(c), cimag(c));
    else
        printf("%lf - i%lf", creal(c), fabs(cimag(c)));
}
double complex* makeW(int n) {
    double complex* w = (double complex*)calloc(n, sizeof(double complex));
    for(int i = 0;i < n;i++) 
        w[i] = cexp(- I * 2.0 * M_PI * (double)i / (double)n);
    return w;
}