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

void QD::qd_mul_d_qd(const QD* a, double b, QD* c) {
    double t[4];
    util::two_prod(a->x[0], b, &(c->x[0]), &t[0]);
    util::two_prod(a->x[1], b, &(c->x[1]), &t[1]);
    util::two_sum(c->x[1], t[0], &(c->x[1]), &t[0]);
    util::two_prod(a->x[2], b, &(c->x[2]), &t[2]);
    util::three_sum(c->x[2], t[1], t[0], &(c->x[2]), &t[0], &t[1]);
    c->x[3] = a->x[3] * b;
    util::three_sum(c->x[3], t[2], t[0], &(c->x[3]), &t[0]);
    t[0] += t[1];
    c->renormalize(t[0]);
}

void QD::qd_mul_qd_qd(const QD* a, const QD* b, QD* c) {
    double t[13];
    util::two_prod(a->x[0], b->x[0], &(c->x[0]), &(c->x[1])); // 0, 1  
    util::two_prod(a->x[0], b->x[1], &t[0], &t[1]); // 1, 2 
    util::two_prod(a->x[1], b->x[0], &t[2], &t[3]); // 1, 2 
    util::three_sum(c->x[1], t[0], t[2], &(c->x[1]), &t[0], &t[2]); // 1, 2, 3 
    util::two_prod(a->x[0], b->x[2], &(c->x[2]), &t[4]); // 2, 3 
    util::two_prod(a->x[1], b->x[1], &t[5], &t[6]); // 2, 3
    util::two_prod(a->x[2], b->x[0], &t[7], &t[8]); // 2, 3 
    util::six_three_sum(t[1], t[3], t[0], c->x[2], t[5], t[7], &(c->x[2]), &t[0], &t[1]); // 2, 3, 4 
    util::two_prod(a->x[0], b->x[3], &(c->x[3]), &t[3]); // 3, 4
    util::two_prod(a->x[1], b->x[2], &t[5], &t[7]);  // 3, 4
    util::two_prod(a->x[2], b->x[1], &t[9], &t[10]);  // 3, 4
    util::two_prod(a->x[3], b->x[0], &t[11], &t[12]);  // 3, 4
    util::nine_two_sum(t[2], t[4], t[6], t[8], t[0], c->x[3], t[5], t[9], t[11], &(c->x[3]), &t[0]); // 3, 4 
    t[0] = a->x[1] * b->x[3] + a->x[2] * b->x[2] + a->x[3] * b->x[1] + t[1] + t[3] + t[7] + t[10] + t[12] + t[0];
    c->renormalize(t[0]);
}