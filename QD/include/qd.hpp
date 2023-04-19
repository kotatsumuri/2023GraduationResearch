#pragma once

class qd {
    public:
        double x[4];
        void renormalize();
        void renormalize(double a);
        static void qd_add_d_qd(qd a, double b, qd* c);
};