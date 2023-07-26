#include "qd.hpp"

namespace QD {
void zero(qd a) {
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
    a[3] = 0.0;
}

void init(qd a, double x0) {
    a[0] = x0;
    a[1] = 0.0;
    a[2] = 0.0;
    a[3] = 0.0;
}

void init(qd a, double x0, double x1, double x2, double x3) {
    a[0] = x0;
    a[1] = x1;
    a[2] = x2;
    a[3] = x3;
}

void copy(const qd a, qd b) {
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
    b[3] = a[3];
}

void minus(const qd a, qd b) {
    b[0] = -a[0];
    b[1] = -a[1];
    b[2] = -a[2];
    b[3] = -a[3];
}

void minus(qd a) {
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
    a[3] = -a[3];
}

}  // namespace QD