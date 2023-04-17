#pragma once

void decimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w);
void invDecimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w);
void cooleyTukey(int n, int m, double complex* x, double complex* w);
void invCooleyTukey(int n, int m, double complex* x, double complex* w);
void stockham(int n, int p, double complex* x, double complex* y, double complex* w);
void invStockham(int n, int p, double complex* x, double complex* y, double complex* w);
void realFFT(int n, int m, double* x, double complex* y, double complex* w);
void invRealFFT(int n, int m, double* x, double complex* y, double complex* w);
void splitRadixFFT(int initN, int n, double complex* x, double complex* w);
void invSplitRadixFFT(int initN, int n, double complex* x, double complex* w);