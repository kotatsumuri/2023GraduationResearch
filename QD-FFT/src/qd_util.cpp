#include <bitset>
#include <random>

#include "qd.hpp"
#include "util_calc.hpp"

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

bool equal(const qd a, const qd b) {
    return !(a[0] != b[0] || a[1] != b[1] || a[2] != b[2] || a[3] != b[3]);
}

void rand(qd a) {
    static std::random_device rnd;
    static std::mt19937 engine{rnd()};
    static std::uniform_real_distribution<> dist{0, 1};
    static const double m_const = 4.6566128730773926e-10;
    double m                    = 1;
    QD::zero(a);
    double d;

    for (int i = 0; i < 7; i++, m *= m_const) {
        d = dist(engine) * m;
        QD::add(a, d, a);
    }
}

int ulp(qd a) {
    int ret = util::calc::ulp(a[3]);
    int i = 2;
    while(ret <= -1024) {
        ret = util::calc::ulp(a[i]);
        i--;
    }
    return ret;
}

void abs(const qd a, qd b) {
    if(a[0] >= 0)
        QD::copy(a, b);
    else
        QD::minus(a, b);    
}

}  // namespace QD