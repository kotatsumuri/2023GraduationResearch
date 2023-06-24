#include "qd.hpp"

#include <cmath>

#include "util_calc.hpp"
using namespace util::calc;

namespace QD {
    inline void zero(qd a) {
        a[0] = 0.0;
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;
    }

    inline void init(qd a, double x0) {
        a[0] = x0;
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;
    }

    inline void init(qd a, double x0, double x1, double x2, double x3) {
        a[0] = x0;
        a[1] = x1;
        a[2] = x2;
        a[3] = x3;
    }

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

    void mul(const qd a, const qd b, qd p) {
        double t[17];
        two_prod(a[0], b[0], p, p + 1);                // 0, 1
        two_prod(a[0], b[1], t, t + 1);                // 1, 2
        two_prod(a[1], b[0], t + 2, t + 3);            // 1, 2
        three_sum(p[1], t[0], t[2], p + 1, t, t + 2);  // 1, 2, 3
        two_prod(a[0], b[2], p + 2, t + 4);            // 2, 3
        two_prod(a[1], b[1], t + 5, t + 6);            // 2, 3
        two_prod(a[2], b[0], t + 7, t + 8);            // 2, 3
        six_three_sum(t[1], t[3], t[0], p[2], t[5], t[7], p + 2, t,
                      t + 1);                          // 2, 3, 4
        two_prod(a[0], b[3], p + 3, t + 3);            // 3, 4
        two_prod(a[1], b[2], t + 5, t + 7);            // 3, 4
        two_prod(a[2], b[1], t + 9, t + 10);           // 3, 4
        two_prod(a[3], b[0], t + 11, t + 12);          // 3, 4
        nine_two_sum(t[2], t[4], t[6], t[8], t[0], p[3], t[5], t[9], t[11],
                     p + 3, t);                        // 3, 4
        t[0] = a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + t[1] + t[3] + t[7] +
               t[10] + t[12] + t[0];
        renormalize(p, t[0]);
    }

    void mul(const qd a, double b, qd p) {
        double t[4];
        two_prod(a[0], b, p, t);
        two_prod(a[1], b, p + 1, t + 1);
        two_sum(p[1], t[0], p + 1, t);
        two_prod(a[2], b, p + 2, t + 2);
        three_sum(p[2], t[1], t[0], p + 2, t, t + 1);
        p[3] = a[3] * b;
        three_sum(p[3], t[2], t[0], p + 3, t);
        t[0] += t[1];
        renormalize(p, t[0]);
    }

    void mul(double a, double b, qd p) {
        two_prod(a, b, p, p + 1);
        p[2] = p[3] = 0;
    }

    void mul_pwr2(const qd a, double b, qd p) {
        p[0] = a[0] * b;
        p[1] = a[1] * b;
        p[2] = a[2] * b;
        p[3] = a[3] * b;
    }

    void div(const qd a, const qd b, qd d) {
        qd t0, t1, t2;
        d[0] = a[0] / b[0];
        mul(b, d[0], t0);
        sub(a, t0, t1);
        d[1] = t1[0] / b[0];
        mul(b, d[1], t0);
        sub(t1, t0, t2);
        d[2] = t2[0] / b[0];
        mul(b, d[2], t0);
        sub(t2, t0, t1);
        d[3] = t1[0] / b[0];
        mul(b, d[3], t0);
        sub(t1, t0, t2);
        t2[0] = t2[0] / b[0];
        renormalize(d, t2[0]);
    }

    void div(const qd a, double b, qd d) {
        qd t0, t1, t2;
        d[0] = a[0] / b;
        mul(b, d[0], t0);
        sub(a, t0, t1);
        d[1] = t1[0] / b;
        mul(b, d[1], t0);
        sub(t1, t0, t2);
        d[2] = t2[0] / b;
        mul(b, d[2], t0);
        sub(t2, t0, t1);
        d[3] = t1[0] / b;
        mul(b, d[3], t0);
        sub(t1, t0, t2);
        t2[0] = t2[0] / b;
        renormalize(d, t2[0]);
    }

    void div(double a, const qd b, qd d) {
        qd t0, t1, t2;
        d[0] = a / b[0];
        mul(b, d[0], t0);
        sub(a, t0, t1);
        d[1] = t1[0] / b[0];
        mul(b, d[1], t0);
        sub(t1, t0, t2);
        d[2] = t2[0] / b[0];
        mul(b, d[2], t0);
        sub(t2, t0, t1);
        d[3] = t1[0] / b[0];
        mul(b, d[3], t0);
        sub(t1, t0, t2);
        t2[0] = t2[0] / b[0];
        renormalize(d, t2[0]);
    }

    void div(double a, double b, qd d) {
        qd t0, t1, t2;
        d[0] = a / b;
        mul(b, d[0], t0);
        sub(a, t0, t1);
        d[1] = t1[0] / b;
        mul(b, d[1], t0);
        sub(t1, t0, t2);
        d[2] = t2[0] / b;
        mul(b, d[2], t0);
        sub(t2, t0, t1);
        d[3] = t1[0] / b;
        mul(b, d[3], t0);
        sub(t1, t0, t2);
        t2[0] = t2[0] / b;
        renormalize(d, t2[0]);
    }

    void sqrt(const qd a, qd b) {
        if (a[0] != 0.0 || a[1] != 0.0 || a[2] != 0.0 || a[3] == 0.0) {
            qd h;
            qd t0, t1;
            mul_pwr2(a, 0.5, h);
            init(b, 1 / std::sqrt(a[0]));

            mul(b, b, t0);
            mul(h, t0, t1);
            sub(0.5, t1, t0);
            mul(b, t0, t1);
            add(b, t1, t0);

            mul(t0, t0, b);
            mul(h, b, t1);
            sub(0.5, t1, b);
            mul(t0, b, t1);
            add(t0, t1, b);

            mul(b, b, t0);
            mul(h, t0, t1);
            sub(0.5, t1, t0);
            mul(b, t0, t1);
            add(b, t1, t0);

            mul(a, t0, b);
        } else {
            zero(b);
        }
    }
}  // namespace QD