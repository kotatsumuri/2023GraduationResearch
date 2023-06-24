#include "qd.hpp"

#include "util_calc.hpp"
using namespace util::calc;

namespace QD {
    void renormalize(qd a) {
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

    void renormalize(qd a, double b) {
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

    void add(const qd a, const qd b, qd s) {
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

    void add(const qd a, double b, qd s) {
        double e;
        two_sum(a[0], b, s, &e);
        two_sum(a[1], e, s + 1, &e);
        two_sum(a[2], e, s + 2, &e);
        two_sum(a[3], e, s + 3, &e);
        renormalize(s, e);
    }

    void sub(const qd a, const qd b, qd s) {
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

    void sub(const qd a, double b, qd s) {
        double e;
        two_sum(a[0], -b, s, &e);
        two_sum(a[1], e, s + 1, &e);
        two_sum(a[2], e, s + 2, &e);
        two_sum(a[3], e, s + 3, &e);
        renormalize(s, e);
    }

    void sub(double a, const qd b, qd s) {
        double e;
        two_sum(-b[0], a, s, &e);
        two_sum(-b[1], e, s + 1, &e);
        two_sum(-b[2], e, s + 2, &e);
        two_sum(-b[3], e, s + 3, &e);
        renormalize(s, e);
    }
}  // namespace QD