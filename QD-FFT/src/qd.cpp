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
}  // namespace QD