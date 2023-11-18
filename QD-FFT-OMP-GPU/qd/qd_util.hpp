#pragma once

#include "qd.hpp"

inline int log_ufp(double a) {
    if (a == 0)
        return -1024;
    return std::max(int(floor(log2(fabs(a)))), -1023);
}
inline int log_ulp(double a) { return log_ufp(a) - 52; }
inline int log_ufp(qd a) {
    int ret = log_ufp(a[0]);
    int i   = 1;
    while (ret <= -1024 || i < 4) ret = log_ufp(a[i++]);
    return ret;
}
inline int log_ulp(const qd a) {
    int ret = log_ulp(a[3]);
    int i   = 2;
    while (ret <= -1024 - 52) ret = log_ulp(a[i--]);
    return ret;
}
