#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void DFT(int n, double complex* x, double complex *y, double complex *w);
void invDFT(int n, double complex* x, double complex *y, double complex *w);