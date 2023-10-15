#pragma once

#include <cmath>

#include "qd.hpp"

inline void sqrt(const qd a, qd b) {
    if (a[0] != 0.0 || a[1] != 0.0 || a[2] != 0.0 || a[3] != 0.0) {
        qd h;
        qd t0, t1;
        mul_pwr2(a, 0.5, h);
        init(b, 1.0 / std::sqrt(a[0]));

        sqr(b, t0);
        mul(h, t0, t1);
        sub(0.5, t1, t0);
        mul(b, t0, t1);
        add(b, t1, t0);

        sqr(t0, b);
        mul(h, b, t1);
        sub(0.5, t1, b);
        mul(t0, b, t1);
        add(t0, t1, b);

        sqr(b, t0);
        mul(h, t0, t1);
        sub(0.5, t1, t0);
        mul(b, t0, t1);
        add(b, t1, t0);

        mul(a, t0, b);
    } else {
        zero(b);
    }
}

inline void sqrt(double a, qd b) {
    if (a != 0) {
        qd h;
        qd t0, t1;
        init(h, a * 0.5);
        init(b, 1 / std::sqrt(a));

        sqr(b, t0);
        mul(h, t0, t1);
        sub(0.5, t1, t0);
        mul(b, t0, t1);
        add(b, t1, t0);

        sqr(t0, b);
        mul(h, b, t1);
        sub(0.5, t1, b);
        mul(t0, b, t1);
        add(t0, t1, b);

        sqr(b, t0);
        mul(h, t0, t1);
        sub(0.5, t1, t0);
        mul(b, t0, t1);
        add(b, t1, t0);

        mul(t0, a, b);
    } else {
        zero(b);
    }
}