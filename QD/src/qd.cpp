#include <iostream>
#include "qd.hpp"
#include "util.hpp"

std::ostream& operator << (std::ostream& os, const QD& a) {
    os << util::to_reg_str(a.x[0]) << "|" << util::to_reg_str(a.x[1]) << "|" << util::to_reg_str(a.x[2]) << "|" << util::to_reg_str(a.x[3]);
    return os;
}

QD::QD() {
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
}
QD::QD(double x0, double x1, double x2, double x3) {
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    x[3] = x3;
    renormalize();
}

void QD::renormalize() {
    double s, e, t[5];
    int k, i;

    t[4] = 0.0;
    util::quick_two_sum(x[2], x[3], &s,    &t[3]);
    util::quick_two_sum(x[1], s,    &s,    &t[2]);
    util::quick_two_sum(x[0], s,    &t[0], &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        util::quick_two_sum(s, t[i], &s, &e);
        if(e != 0) {
            x[k++] = s;
            s = e;
        }
    }

    if(k < 4) {
        x[k++] = s;
        while(k < 4)
            x[k++] = 0.0;
    }
}

void QD::renormalize(double a) {
    double s, e, t[5];
    int k, i;

    util::quick_two_sum(x[3], a, &s,    &t[4]);
    util::quick_two_sum(x[2], s, &s,    &t[3]);
    util::quick_two_sum(x[1], s, &s,    &t[2]);
    util::quick_two_sum(x[0], s, &t[0], &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        util::quick_two_sum(s, t[i], &s, &e);
        if(e != 0) {
            x[k++] = s;
            s = e;
        }
    }

    if(k < 4) {
        x[k++] = s;
        while(k < 4)
            x[k++] = 0.0;
    }
}

void QD::qd_add_d_qd(const QD* a, double b, QD* c) {
    double e;
    
    util::two_sum(a->x[0], b, &(c->x[0]), &e);
    util::two_sum(a->x[1], e, &(c->x[1]), &e);
    util::two_sum(a->x[2], e, &(c->x[2]), &e);
    util::two_sum(a->x[3], e, &(c->x[3]), &e);

    c->renormalize(e);
}

void QD::qd_add_qd_qd(const QD* a, const QD* b, QD* c) {
    double t[4];
    util::two_sum(a->x[0],  b->x[0], &(c->x[0]), &t[0]);
    util::two_sum(a->x[1],  b->x[1], &(c->x[1]), &t[1]);
    util::two_sum(c->x[1], t[0],   &(c->x[1]), &t[0]);
    util::two_sum(a->x[2],  b->x[2], &(c->x[2]), &t[2]);
    util::three_sum(c->x[2], t[1], t[0], &(c->x[2]), &t[0], &t[1]);
    util::two_sum(a->x[3],  b->x[3], &(c->x[3]), &t[3]);
    util::three_sum(c->x[3], t[2], t[0], &(c->x[3]), &t[0]);
    t[1] += t[0] + t[3];
    c->renormalize(t[1]);
}