#pragma once
#include "qd.hpp"
#include "qd_init.hpp"

inline void renormalize(qd a) {
    double s, e, t[5];
    int k, i;

    t[4] = 0.0;
    quick_two_sum(a[2], a[3], &s, t + 3);
    quick_two_sum(a[1], s, &s, t + 2);
    quick_two_sum(a[0], s, t, t + 1);

    s = t[0];
    k = 0;
    for (i = 1; i <= 4; i++) {
        quick_two_sum(s, t[i], &s, &e);
        if (e != 0) {
            a[k++] = s;
            s      = e;
        }
    }

    if (k < 4) {
        a[k++] = s;
        while (k < 4) a[k++] = 0.0;
    }
}

inline void renormalize(qd restrict a, double b) {
    double s, e, t[5];
    int k, i;

    quick_two_sum(a[3], b, &s, t + 4);
    quick_two_sum(a[2], s, &s, t + 3);
    quick_two_sum(a[1], s, &s, t + 2);
    quick_two_sum(a[0], s, t, t + 1);

    s = t[0];
    k = 0;
    for (i = 1; i <= 4; i++) {
        quick_two_sum(s, t[i], &s, &e);
        if (e != 0) {
            a[k++] = s;
            s      = e;
        }
    }

    if (k < 4) {
        a[k++] = s;
        while (k < 4) a[k++] = 0.0;
    }
}

inline void add(const qd a, const qd b, qd restrict s) {
    double t[4];
    two_sum(a[0], b[0], s, t);
    two_sum(a[1], b[1], s + 1, t + 1);
    two_sum(s[1], t[0], s + 1, t);
    two_sum(a[2], b[2], s + 2, t + 2);
    three_sum(s[2], t[1], t[0], s + 2, t, t + 1);
    two_sum(a[3], b[3], s + 3, t + 3);
    three_sum(s[3], t[2], t[0], s + 3, t);
    t[1] += t[0] + t[3];
    renormalize(s, t[1]);
}

inline void add(const qd a, double b, qd restrict s) {
    double e;
    two_sum(a[0], b, s, &e);
    two_sum(a[1], e, s + 1, &e);
    two_sum(a[2], e, s + 2, &e);
    two_sum(a[3], e, s + 3, &e);
    renormalize(s, e);
}

inline void sub(const qd a, const qd b, qd restrict s) {
    double t[4];
    two_sum(a[0], -b[0], s, t);
    two_sum(a[1], -b[1], s + 1, t + 1);
    two_sum(s[1], t[0], s + 1, t);
    two_sum(a[2], -b[2], s + 2, t + 2);
    three_sum(s[2], t[1], t[0], s + 2, t, t + 1);
    two_sum(a[3], -b[3], s + 3, t + 3);
    three_sum(s[3], t[2], t[0], s + 3, t);
    t[1] += t[0] + t[3];
    renormalize(s, t[1]);
}

inline void sub(const qd a, double b, qd restrict s) {
    double e;
    two_sum(a[0], -b, s, &e);
    two_sum(a[1], e, s + 1, &e);
    two_sum(a[2], e, s + 2, &e);
    two_sum(a[3], e, s + 3, &e);
    renormalize(s, e);
}

inline void sub(double a, const qd b, qd restrict s) {
    double e;
    two_sum(-b[0], a, s, &e);
    two_sum(-b[1], e, s + 1, &e);
    two_sum(-b[2], e, s + 2, &e);
    two_sum(-b[3], e, s + 3, &e);
    renormalize(s, e);
}

inline void mul(const qd a, const qd b, qd restrict p) {
    double t[17];
    two_prod(a[0], b[0], p, p + 1);                // 0, 1
    two_prod(a[0], b[1], t, t + 1);                // 1, 2
    two_prod(a[1], b[0], t + 2, t + 3);            // 1, 2
    three_sum(p[1], t[0], t[2], p + 1, t, t + 2);  // 1, 2, 3
    two_prod(a[0], b[2], p + 2, t + 4);            // 2, 3
    two_prod(a[1], b[1], t + 5, t + 6);            // 2, 3
    two_prod(a[2], b[0], t + 7, t + 8);            // 2, 3
    six_three_sum(t[1], t[3], t[0], p[2], t[5], t[7], p + 2, t,
                  t + 1);  // 2, 3, 4
#ifdef SLOPPY_MUL
    p[3] = a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0] + t[0] + t[2] +
           t[4] + t[6] + t[8];
    renormalize(p, t[1]);
#else
    two_prod(a[0], b[3], p + 3, t + 3);    // 3, 4
    two_prod(a[1], b[2], t + 5, t + 7);    // 3, 4
    two_prod(a[2], b[1], t + 9, t + 10);   // 3, 4
    two_prod(a[3], b[0], t + 11, t + 12);  // 3, 4
    nine_two_sum(t[2], t[4], t[6], t[8], t[0], p[3], t[5], t[9], t[11], p + 3,
                 t);  // 3, 4
    t[0] = a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + t[1] + t[3] + t[7] +
           t[10] + t[12] + t[0];
    renormalize(p, t[0]);
#endif
}

void mul(const qd a, double b, qd restrict p) {
    double t[4];
    two_prod(a[0], b, p, t);
    two_prod(a[1], b, p + 1, t + 1);
    two_sum(p[1], t[0], p + 1, t);
    two_prod(a[2], b, p + 2, t + 2);
    three_sum(p[2], t[1], t[0], p + 2, t, t + 1);
    p[3] = a[3] * b;
    three_sum(p[3], t[2], t[0], p + 3, t);
    t[0] += t[1];
    renormalize(p, t[0]);
}

inline void mul_pwr2(const qd a, double b, qd restrict p) {
    p[0] = a[0] * b;
    p[1] = a[1] * b;
    p[2] = a[2] * b;
    p[3] = a[3] * b;
}

inline void div_pwr2(const qd a, double b, qd restrict p) {
    p[0] = a[0] / b;
    p[1] = a[1] / b;
    p[2] = a[2] / b;
    p[3] = a[3] / b;
}

inline void sqr(const qd a, qd restrict p) {
    double t[17];
    two_prod(a[0], a[0], p, p + 1);                // 0, 1
    two_prod(a[0], a[1], t, t + 1);                // 1, 2
    two_prod(a[1], a[0], t + 2, t + 3);            // 1, 2
    three_sum(p[1], t[0], t[2], p + 1, t, t + 2);  // 1, 2, 3
    two_prod(a[0], a[2], p + 2, t + 4);            // 2, 3
    two_prod(a[1], a[1], t + 5, t + 6);            // 2, 3
    two_prod(a[2], a[0], t + 7, t + 8);            // 2, 3
    six_three_sum(t[1], t[3], t[0], p[2], t[5], t[7], p + 2, t,
                  t + 1);                  // 2, 3, 4
    two_prod(a[0], a[3], p + 3, t + 3);    // 3, 4
    two_prod(a[1], a[2], t + 5, t + 7);    // 3, 4
    two_prod(a[2], a[1], t + 9, t + 10);   // 3, 4
    two_prod(a[3], a[0], t + 11, t + 12);  // 3, 4
    nine_two_sum(t[2], t[4], t[6], t[8], t[0], p[3], t[5], t[9], t[11], p + 3,
                 t);  // 3, 4
    t[0] = a[1] * a[3] + a[2] * a[2] + a[3] * a[1] + t[1] + t[3] + t[7] +
           t[10] + t[12] + t[0];
    renormalize(p, t[0]);
}

inline void fabs(const qd a, qd restrict b) {
    if (a[0] >= 0)
        copy(a, b);
    else
        minus(a, b);
}