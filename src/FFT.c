#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "FFT.h"
#include "Utility.h"

void decimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w) {
    int i;
    if(n <= 1)
        return;

    for(i = 0;i < n / 2;i++) {
        y[i] = x[i] + x[i + n / 2];
        y[i + n / 2] = (x[i] - x[i + n / 2]) * w[i * initN / n];        
    }

    decimationInFrequency(initN, n / 2, y, x, w);
    decimationInFrequency(initN, n / 2, &y[n / 2], x, w);

    for(i = 0;i < n / 2;i++) {
        x[2 * i] = y[i];
        x[2 * i + 1] = y[i + n / 2];
    }
}

void invDecimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w) {
    int i;
    if(n <= 1) {
        x[0] /= initN;
        return;
    }

    for(i = 0;i < n / 2;i++) {
        y[i] = x[i] + x[i + n / 2];
        y[i + n / 2] = (x[i] - x[i + n / 2]) / w[i * initN / n];        
    }

    invDecimationInFrequency(initN, n / 2, y, x, w);
    invDecimationInFrequency(initN, n / 2, &y[n / 2], x, w);

    for(i = 0;i < n / 2;i++) {
        x[2 * i] = y[i];
        x[2 * i + 1] = y[i + n / 2];
    }
}

void cooleyTukey(int n, int m, double complex* x, double complex* w) {
    int i, j, k, l;
    double complex temp, W;
    l = n;
    for(k = 1;k <= m;k++) {
        l /= 2;
        for(j = 1;j <= l;j++) {
            W = w[(j - 1) * n / (l * 2)];
            for(i = j;i <= n;i += l * 2) {
                temp = x[i - 1] - x[i + l - 1];
                x[i - 1] = x[i - 1] + x[i + l - 1];
                x[i + l - 1] = temp * W;
            }
        }
    }

    bitReversed(n, x);
}

void invCooleyTukey(int n, int m, double complex* x, double complex* w) {
    int i, j, k, l;
    double complex temp, W;
    l = n;
    for(k = 1;k <= m;k++) {
        l /= 2;
        for(j = 1;j <= l;j++) {
            W = w[(j - 1) * n / (l * 2)];
            for(i = j;i <= n;i += l * 2) {
                temp = x[i - 1] - x[i + l - 1];
                x[i - 1] = x[i - 1] + x[i + l - 1];
                x[i + l - 1] = temp / W;
            }
        }
    }

    for(i = 0;i < n;i++)
        x[i] /= n;

    bitReversed(n, x);
}

void stockham(int n, int p, double complex* x, double complex* y, double complex* w) {
    int l = n / 2, m = 1, t, j, k, flag = 0;
    double complex c0, c1;
    double complex** X = (double complex**)calloc(2, sizeof(double complex*));
    X[0] = (double complex*)calloc(n, sizeof(double complex));
    X[1] = (double complex*)calloc(n, sizeof(double complex));
    for(j = 0;j < n;j++)
        X[0][j] = x[j];

    for(t = 1;t <= p;t++) {
        for(j = 0;j <= l - 1;j++) {
            for(k = 0;k <= m - 1;k++) {
                c0 = X[flag][k + j * m];
                c1 = X[flag][k + j * m + l * m];
                X[!flag][k + 2 * j * m] = c0 + c1;
                X[!flag][k + 2 * j * m + m] = w[j * n / (2 * l)] * (c0 - c1);
            }
        }
        l /= 2;
        m *= 2;
        flag = !flag;
    }

    for(j = 0;j < n;j++)
        y[j] = X[flag][j];
}

void invStockham(int n, int p, double complex* x, double complex* y, double complex* w) {
    int l = n / 2, m = 1, t, j, k, flag = 0;
    double complex c0, c1;
    double complex** X = (double complex**)calloc(2, sizeof(double complex*));
    X[0] = (double complex*)calloc(n, sizeof(double complex));
    X[1] = (double complex*)calloc(n, sizeof(double complex));
    for(j = 0;j < n;j++)
        X[0][j] = x[j];

    for(t = 1;t <= p;t++) {
        for(j = 0;j <= l - 1;j++) {
            for(k = 0;k <= m - 1;k++) {
                c0 = X[flag][k + j * m];
                c1 = X[flag][k + j * m + l * m];
                X[!flag][k + 2 * j * m] = c0 + c1;
                X[!flag][k + 2 * j * m + m] = (c0 - c1) / w[j * n / (2 * l)];
            }
        }
        l /= 2;
        m *= 2;
        flag = !flag;
    }

    for(j = 0;j < n;j++)
        y[j] = X[flag][j] / n;
}

void realFFT(int n, int m, double* x, double complex* y, double complex* w) {
    int i;
    double complex* a = (double complex*)calloc(n / 2, sizeof(double complex));
    for(i = 0;i < n / 2;i++)
        a[i] = x[2 * i] + I * x[2 * i + 1];
    double complex *w2 = (double complex*)calloc(n / 2, sizeof(double complex));
    for(i = 0;i < n / 2;i++)
        w2[i] = w[2 * i];

    cooleyTukey(n / 2, m, a, w2);
    
    y[0] = creal(a[0]) + cimag(a[0]);
    y[n / 4] = a[n / 4];
    y[n / 2] = creal(a[0]) - cimag(a[0]);
    
    double complex temp;
    for(i = 1;i < n / 4;i++) {
        temp = (a[i] - conj(a[n / 2 - i])) * (1 + I * w[i]) * 0.5;
        y[i] = a[i] - temp;
        y[n - i] = conj(y[i]);
        y[n / 2 + i] = conj(a[n / 2 - i]) + temp;
        y[n / 2 - i] = conj(y[n / 2 + i]);
    }
}

void invRealFFT(int n, int m, double* x, double complex* y, double complex* w) {
    int i;
    double complex* a = (double complex*)calloc(n / 2, sizeof(double complex));
    for(i = 0;i < n / 2;i++)
        a[i] = x[2 * i] + I * x[2 * i + 1];
    double complex *w2 = (double complex*)calloc(n / 2, sizeof(double complex));
    for(i = 0;i < n / 2;i++)
        w2[i] = w[2 * i];

    invCooleyTukey(n / 2, m, a, w2);
    
    y[0] = (creal(a[0]) + cimag(a[0])) / 2;
    y[n / 4] = a[n / 4] / 2;
    y[n / 2] = (creal(a[0]) - cimag(a[0])) / 2;
    
    double complex temp;
    for(i = 1;i < n / 4;i++) {
        temp = (a[i] - conj(a[n / 2 - i])) * (1 + I / w[i]) * 0.5;
        y[i] = (a[i] - temp) / 2;
        y[n - i] = conj(y[i]);
        y[n / 2 + i] = (conj(a[n / 2 - i]) + temp) / 2;
        y[n / 2 - i] = conj(y[n / 2 + i]);
    }
}