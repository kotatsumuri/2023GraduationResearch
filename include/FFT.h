#pragma once

void decimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w);
void invDecimationInFrequency(int initN, int n, double complex* x, double complex* y, double complex* w);
void cooleyTukey(int n, int m, double complex* x, double complex* w);
void invCooleyTukey(int n, int m, double complex* x, double complex* w);