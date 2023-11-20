#pragma once
#include "../qd/qd.hpp"

inline void butterfly(qd x0, qd ix0, qd x1, qd ix1, qd y0, qd iy0, qd y1, qd iy1, const qd a, const qd b) {
    add(x0, x1, y0);
    add(ix0, ix1, iy0);

    sub(x0, x1, y1);
    sub(ix0, ix1, iy1);

    mul(y1, a, x0);
    mul(y1, b, ix0);
    mul(iy1, a, x1);
    mul(iy1, b, ix1);

    add(x0, ix1, y1);
    sub(x1, ix0, iy1);
}

inline void inv_butterfly(qd x0, qd ix0, qd x1, qd ix1, qd y0, qd iy0, qd y1, qd iy1, const qd a, const qd b) {
    add(x0, x1, y0);
    add(ix0, ix1, iy0);

    sub(x0, x1, y1);
    sub(ix0, ix1, iy1);

    mul(y1, a, x0);
    mul(y1, b, ix0);
    mul(iy1, a, x1);
    mul(iy1, b, ix1);

    sub(x0, ix1, y1);
    add(x1, ix0, iy1);
}