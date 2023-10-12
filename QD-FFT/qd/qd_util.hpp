#pragma once

#include "qd.hpp"

inline int log_ufp(double a) {
    return std::max(int(floor(log2(fabs(a)))), -1023);
}
inline int log_ulp(double a) { return log_ufp(a) - 52; }
inline int ulp(qd a) {
    int ret = log_ulp(a[3]);
    int i   = 2;
    while (ret <= -1024) ret = log_ulp(a[i--]);
    return ret;
}