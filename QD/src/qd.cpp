#include "qd.hpp"
#include "utility.hpp"

void qd::renormalize() {
    double s, e, t[5];
    int k, i;

    t[4] = 0.0;
    utility::quick_two_sum(x[2], x[3], &s,    &t[3]);
    utility::quick_two_sum(x[1], s,    &s,    &t[2]);
    utility::quick_two_sum(x[0], s,    &t[0], &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        utility::quick_two_sum(s, t[i], &s, &e);
        if(e != 0) {
            x[k] = s;
            s = e;
            k++;
        }
    }
}

void qd::renormalize(double a) {
    double s, e, t[5];
    int k, i;

    utility::quick_two_sum(x[3], a, &s,    &t[4]);
    utility::quick_two_sum(x[2], s, &s,    &t[3]);
    utility::quick_two_sum(x[1], s, &s,    &t[2]);
    utility::quick_two_sum(x[0], s, &t[0], &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        utility::quick_two_sum(s, t[i], &s, &e);
        if(e != 0) {
            x[k] = s;
            s = e;
            k++;
        }
    }
}

void qd::qd_add_d_qd(qd a, double b, qd* c) {
    double e;
    
    utility::two_sum(a.x[0], b, &(c->x[0]), &e);
    utility::two_sum(a.x[1], e, &(c->x[1]), &e);
    utility::two_sum(a.x[2], e, &(c->x[2]), &e);
    utility::two_sum(a.x[3], e, &(c->x[3]), &e);

    c->renormalize(e);
}

void qd::qd_add_qd_qd(qd a, qd b, qd* c) {
    double t[4];
    utility::two_sum(a.x[0],  b.x[0], &(c->x[0]), &t[0]);
    utility::two_sum(a.x[1],  b.x[1], &(c->x[1]), &t[1]);
    utility::two_sum(c->x[1], t[0],   &(c->x[1]), &t[0]);
    utility::two_sum(a.x[2],  b.x[2], &(c->x[2]), &t[2]);
    utility::three_sum(c->x[2], t[1], t[0], &(c->x[2]), &t[0], &t[1]);
    utility::two_sum(a.x[3],  b.x[3], &(c->x[3]), &t[3]);
    utility::three_sum(c->x[3], t[2], t[0], &(c->x[3]), &t[0]);
    t[1] += t[0] + t[3];
    c->renormalize(t[1]);
}