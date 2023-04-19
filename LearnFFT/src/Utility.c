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

void bitReversed(int n, double complex* x) {
    int i, j, k;
    double complex temp;
    j = 1;
    for(i = 1;i <= n - 1;i++) {
        if(i < j) {
            temp = x[i - 1];
            x[i - 1] = x[j - 1];
            x[j - 1] = temp;
        }
        k = n / 2;
        while(k < j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}