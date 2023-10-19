#pragma once
#include "../qd/qd.hpp"

void swap(double **a, double **b) {
    double *tmp = *a;
    *a          = *b;
    *b          = tmp;
}

void swap(qd **a, qd **b) {
    qd *tmp = *a;
    *a      = *b;
    *b      = tmp;
}