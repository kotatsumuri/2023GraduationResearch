#pragma once
#include <iostream>

class QD {
    friend std::ostream& operator << (std::ostream& os, const QD& a);
    
    public:
        double x[4];
        QD();
        QD(double x0, double x1, double x2, double x3);
        void renormalize();
        void renormalize(double a);
        static void qd_add_d_qd(const QD* a, double b, QD* c);
        static void qd_add_qd_qd(const QD* a, const QD* b, QD* c);  
        static void qd_prod_d_qd(const QD* a, double b, QD* c);
        static void qd_prod_qd_qd(const QD* a, const QD* b, QD* c);
};