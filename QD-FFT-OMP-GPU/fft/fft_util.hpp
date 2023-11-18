#pragma once
#include "../qd/qd.hpp"

inline void swap(double **a, double **b) {
    double *tmp = *a;
    *a          = *b;
    *b          = tmp;
}

inline void swap(double ***a, double ***b) {
    double **tmp = *a;
    *a           = *b;
    *b           = tmp;
}

inline void swap(qd **a, qd **b) {
    qd *tmp = *a;
    *a      = *b;
    *b      = tmp;
}

inline void swap(qd ***a, qd ***b) {
    qd **tmp = *a;
    *a       = *b;
    *b       = tmp;
}

inline void swap(qd a, qd b) {
    double tmp = a[0];
    a[0]       = b[0];
    b[0]       = tmp;
    tmp        = a[1];
    a[1]       = b[1];
    b[1]       = tmp;
    tmp        = a[2];
    a[2]       = b[2];
    b[2]       = tmp;
    tmp        = a[3];
    a[3]       = b[3];
    b[3]       = tmp;
}