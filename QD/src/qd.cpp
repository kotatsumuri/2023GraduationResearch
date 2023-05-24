#include <iostream>
#include <cmath>
#include "qd.hpp"
#include "util.hpp"

using namespace QD_Lib;

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

QD::QD(double x0) {
    x[0] = x0;
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
    s = util::quick_two_sum(x[2], x[3], &t[3]);
    s = util::quick_two_sum(x[1], s, &t[2]);
    t[0] = util::quick_two_sum(x[0], s, &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        s = util::quick_two_sum(s, t[i], &e);
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

    s = util::quick_two_sum(x[3], a, &t[4]);
    s = util::quick_two_sum(x[2], s, &t[3]);
    s = util::quick_two_sum(x[1], s, &t[2]);
    t[0] = util::quick_two_sum(x[0], s, &t[1]);

    s = t[0];
    k = 0;
    for(i = 1;i <= 4;i++) {
        s = util::quick_two_sum(s, t[i], &e);
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

QD QD::qd_add_d_qd(const QD& a, double b) {
    double e;
    QD c;
    c.x[0] = util::two_sum(a.x[0], b, &e);
    c.x[1] = util::two_sum(a.x[1], e, &e);
    c.x[2] = util::two_sum(a.x[2], e, &e);
    c.x[3] = util::two_sum(a.x[3], e, &e);

    c.renormalize(e);
    return c;
}

QD QD::qd_add_qd_qd(const QD& a, const QD& b) {
    double t[4];
    QD c;
    c.x[0] = util::two_sum(a.x[0], b.x[0], &t[0]);
    c.x[1] = util::two_sum(a.x[1], b.x[1], &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_sum(a.x[2],  b.x[2], &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = util::two_sum(a.x[3],  b.x[3], &t[3]);
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[1] += t[0] + t[3];
    c.renormalize(t[1]);
    return c;
}

QD QD::qd_sub_d_qd(const QD& a, double b) {
    return qd_add_d_qd(a, -b);
}

QD QD::qd_sub_qd_qd(const QD& a, const QD& b) {
    double t[4];
    QD c;
    c.x[0] = util::two_sum(a.x[0], -b.x[0], &t[0]);
    c.x[1] = util::two_sum(a.x[1], -b.x[1], &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_sum(a.x[2],  -b.x[2], &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = util::two_sum(a.x[3],  -b.x[3], &t[3]);
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[1] += t[0] + t[3];
    c.renormalize(t[1]);
    return c;
}

QD QD::d_sub_qd_qd(double a, const QD& b) {
    double t[4];
    QD c;
    c.x[0] = util::two_sum(a, -b.x[0], &t[0]);
    c.x[1] = util::two_sum(-b.x[1], t[0], &t[0]);
    c.x[2] = util::two_sum(-b.x[2], t[0], &t[0]);
    c.x[3] = util::two_sum(-b.x[3], t[0], &t[0]);
    c.renormalize(t[0]);
    return c;
}

QD QD::qd_mul_d_qd(const QD& a, double b) {
    double t[4];
    QD c;
    c.x[0] = util::two_prod(a.x[0], b, &t[0]);
    c.x[1] = util::two_prod(a.x[1], b, &t[1]);
    c.x[1] = util::two_sum(c.x[1], t[0], &t[0]);
    c.x[2] = util::two_prod(a.x[2], b, &t[2]);
    c.x[2] = util::three_sum(c.x[2], t[1], t[0], &t[0], &t[1]);
    c.x[3] = a.x[3] * b;
    c.x[3] = util::three_sum(c.x[3], t[2], t[0], &t[0]);
    t[0] += t[1];
    c.renormalize(t[0]);
    return c;
}

QD QD::qd_mul_qd_qd(const QD& a, const QD& b) {
    double t[17];
    QD c;
    c.x[0] = util::two_prod(a.x[0], b.x[0], &c.x[1]); // 0, 1  
    t[0] = util::two_prod(a.x[0], b.x[1], &t[1]); // 1, 2 
    t[2] = util::two_prod(a.x[1], b.x[0], &t[3]); // 1, 2 
    c.x[1] = util::three_sum(c.x[1], t[0], t[2], &t[0], &t[2]); // 1, 2, 3 
    c.x[2] = util::two_prod(a.x[0], b.x[2], &t[4]); // 2, 3 
    t[5] = util::two_prod(a.x[1], b.x[1], &t[6]); // 2, 3
    t[7] = util::two_prod(a.x[2], b.x[0], &t[8]); // 2, 3 
    c.x[2] = util::six_three_sum(t[1], t[3], t[0], c.x[2], t[5], t[7], &t[0], &t[1]); // 2, 3, 4 
    c.x[3] = util::two_prod(a.x[0], b.x[3], &t[3]); // 3, 4
    t[5] = util::two_prod(a.x[1], b.x[2], &t[7]);  // 3, 4
    t[9] = util::two_prod(a.x[2], b.x[1], &t[10]);  // 3, 4
    t[11] = util::two_prod(a.x[3], b.x[0], &t[12]);  // 3, 4
    c.x[3] = util::nine_two_sum(t[2], t[4], t[6], t[8], t[0], c.x[3], t[5], t[9], t[11], &t[0]); // 3, 4 
    t[0] = a.x[1] * b.x[3] + a.x[2] * b.x[2] + a.x[3] * b.x[1] + t[1] + t[3] + t[7] + t[10] + t[12] + t[0];
    c.renormalize(t[0]);
    return c;
}

QD QD::qd_mul_pwr2_qd(const QD& a, double b) {
    return QD(a.x[0] * b, a.x[1] * b, a.x[2] * b, a.x[3] * b);
}

QD QD::d_mul_d_qd(double a, double b) {
    QD c;
    c.x[0] = util::two_prod(a, b, &c.x[1]);
    c.renormalize();
    return c;
}

QD QD::qd_div_d_qd(const QD& a, double b) {
    double q;
    QD t0, t1;
    QD c;
    c.x[0] = a.x[0] / b;
    t0 = d_mul_d_qd(b, c.x[0]);
    t0 = qd_sub_qd_qd(a, t0);
    c.x[1] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[1]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[2] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[2]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[3] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[3]);
    t0 = qd_sub_qd_qd(t0, t1);
    q = t0.x[0] / b;
    c.renormalize(q);
    return c;
}

QD QD::qd_div_qd_qd(const QD& a, const QD& b) {
    double q;
    QD t0, t1;
    QD c;
    c.x[0] = a.x[0] / b.x[0];
    t0 = qd_mul_d_qd(b, c.x[0]);
    t0 = qd_sub_qd_qd(a, t0);
    c.x[1] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[1]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[2] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[2]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[3] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[3]);
    t0 = qd_sub_qd_qd(t0, t1);
    q = t0.x[0] / b.x[0];
    c.renormalize(q);
    return c;
}

QD QD::d_div_qd_qd(double a, const QD& b) {
    double q;
    QD t0, t1;
    QD c;
    c.x[0] = a / b.x[0];
    t0 = qd_mul_d_qd(b, c.x[0]);
    t0 = d_sub_qd_qd(a, t0);
    c.x[1] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[1]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[2] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[2]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[3] = t0.x[0] / b.x[0];
    t1 = qd_mul_d_qd(b, c.x[3]);
    t0 = qd_sub_qd_qd(t0, t1);
    q = t0.x[0] / b.x[0];
    c.renormalize(q);
    return c;
}

QD QD::d_div_d_qd(double a, double b) {
    double q;
    QD t0, t1;
    QD c;
    c.x[0] = a / b;
    t0 = d_mul_d_qd(b, c.x[0]);
    t0 = d_sub_qd_qd(a, t0);
    c.x[1] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[1]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[2] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[2]);
    t0 = qd_sub_qd_qd(t0, t1);
    c.x[3] = t0.x[0] / b;
    t1 = d_mul_d_qd(b, c.x[3]);
    t0 = qd_sub_qd_qd(t0, t1);
    q = t0.x[0] / b;
    c.renormalize(q);
    return c;
}

QD QD::minus(const QD& a) {
    return QD(-a.x[0], -a.x[1], -a.x[2], -a.x[3]);
}

QD QD::pow(const QD& a, int n) {
    QD c = a;
    QD b = QD(1.0, 0, 0, 0);
    while(n != 0) {
        if(n & 1)
            b = qd_mul_qd_qd(b, c);
        c = qd_mul_qd_qd(c, c);
        n >>= 1;
    }
    return b;
}

QD QD::sqrt(const QD& a) {
    if(a.x[0] == 0.0 && a.x[1] == 0.0 && a.x[2] == 0.0 && a.x[3] == 0.0)
        return QD(0.0);

    QD b = QD(1.0 / std::sqrt(a.x[0]), 0, 0, 0);
    QD h = qd_mul_pwr2_qd(a, 0.5);

    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, d_sub_qd_qd(0.5, qd_mul_qd_qd(h, qd_mul_qd_qd(b, b)))));
    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, d_sub_qd_qd(0.5, qd_mul_qd_qd(h, qd_mul_qd_qd(b, b)))));
    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, d_sub_qd_qd(0.5, qd_mul_qd_qd(h, qd_mul_qd_qd(b, b)))));

    return qd_mul_qd_qd(a, b);
}

QD QD::root(const QD& a, int n) {
    QD b = QD(std::pow(a.x[0],-1.0 / n), 0, 0, 0);
    QD c = d_div_d_qd(1.0, (double)n);
    QD h = qd_div_d_qd(a, (double)n);

    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, qd_sub_qd_qd(c, qd_mul_qd_qd(h, pow(b, n)))));
    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, qd_sub_qd_qd(c, qd_mul_qd_qd(h, pow(b, n)))));
    b = qd_add_qd_qd(b, qd_mul_qd_qd(b, qd_sub_qd_qd(c, qd_mul_qd_qd(h, pow(b, n)))));

    return d_div_qd_qd(1.0, b);
}

static const QD _eps = QD(std::ldexp(1.0, -209));
static const QD _log2 = QD(6.931471805599452862e-01,
                    2.319046813846299558e-17,
                    5.707708438416212066e-34,
                    -3.582432210601811423e-50);
static const int n_inv_fact = 15;
static const QD inv_fact[n_inv_fact] = {
  QD( 1.66666666666666657e-01,  9.25185853854297066e-18,
           5.13581318503262866e-34,  2.85094902409834186e-50),
  QD( 4.16666666666666644e-02,  2.31296463463574266e-18,
           1.28395329625815716e-34,  7.12737256024585466e-51),
  QD( 8.33333333333333322e-03,  1.15648231731787138e-19,
           1.60494162032269652e-36,  2.22730392507682967e-53),
  QD( 1.38888888888888894e-03, -5.30054395437357706e-20,
          -1.73868675534958776e-36, -1.63335621172300840e-52),
  QD( 1.98412698412698413e-04,  1.72095582934207053e-22,
           1.49269123913941271e-40,  1.29470326746002471e-58),
  QD( 2.48015873015873016e-05,  2.15119478667758816e-23,
           1.86586404892426588e-41,  1.61837908432503088e-59),
  QD( 2.75573192239858925e-06, -1.85839327404647208e-22,
           8.49175460488199287e-39, -5.72661640789429621e-55),
  QD( 2.75573192239858883e-07,  2.37677146222502973e-23,
          -3.26318890334088294e-40,  1.61435111860404415e-56),
  QD( 2.50521083854417202e-08, -1.44881407093591197e-24,
           2.04267351467144546e-41, -8.49632672007163175e-58),
  QD( 2.08767569878681002e-09, -1.20734505911325997e-25,
           1.70222792889287100e-42,  1.41609532150396700e-58),
  QD( 1.60590438368216133e-10,  1.25852945887520981e-26,
          -5.31334602762985031e-43,  3.54021472597605528e-59),
  QD( 1.14707455977297245e-11,  2.06555127528307454e-28,
           6.88907923246664603e-45,  5.72920002655109095e-61),
  QD( 7.64716373181981641e-13,  7.03872877733453001e-30,
          -7.82753927716258345e-48,  1.92138649443790242e-64),
  QD( 4.77947733238738525e-14,  4.39920548583408126e-31,
          -4.89221204822661465e-49,  1.20086655902368901e-65),
  QD( 2.81145725434552060e-15,  1.65088427308614326e-31,
          -2.87777179307447918e-50,  4.27110689256293549e-67)
};

QD QD::exp(const QD& a) {
    double k = ldexp(1.0, 16);
    QD inv_k = d_div_d_qd(1.0, k);
    int m = std::floor(a.x[0] / _log2.x[0] + 0.5);
    QD r = qd_mul_qd_qd(qd_sub_qd_qd(a, qd_mul_d_qd(_log2, m)), inv_k);
    double thresh = to_double(qd_mul_qd_qd(inv_k, _eps));

    QD t;
    QD p = qd_mul_qd_qd(r, r);
    QD s = qd_add_qd_qd(r, qd_mul_pwr2_qd(p, 0.5));
    int i = 0;

    do {
        p =  qd_mul_qd_qd(p, r);
        t = qd_mul_qd_qd(p, inv_fact[i++]);
        s = qd_add_qd_qd(s, t);
    } while (std::abs(to_double(t)) > thresh && i < 9);

    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_qd_qd(qd_mul_pwr2_qd(s, 2.0), qd_mul_qd_qd(s, s));
    s = qd_add_d_qd(s, 1.0);

    return qd_mul_d_qd(s, ldexp(1, m));
}

QD QD::log(const QD& a) {
    QD b(std::log(a.x[0]));
    b = qd_add_d_qd(qd_add_qd_qd(b, qd_mul_qd_qd(a, exp(-b))), -1.0);
    b = qd_add_d_qd(qd_add_qd_qd(b, qd_mul_qd_qd(a, exp(-b))), -1.0);
    b = qd_add_d_qd(qd_add_qd_qd(b, qd_mul_qd_qd(a, exp(-b))), -1.0);
    return b;
}

/**
 * @brief Calculate cos(k * 2pi / n) n is power of 2.
 */
QD QD::cos(int k, int n) {
    std::cout << k << " " << n << std::endl;
    if(k == 0)
        return QD(1.0);
    if(n == 1)
        return QD(1.0);
    if(n == 2)
        return QD(-1.0);

    int harf_n = n >> 1;

    if(k > harf_n)
        k = n - k;
    
    if(k > (harf_n >> 1))
        return -qd_mul_pwr2_qd(sqrt(qd_add_d_qd(qd_mul_pwr2_qd(cos(harf_n - k, harf_n), 2.0), 2.0)), 0.5);
    else 
        return qd_mul_pwr2_qd(sqrt(qd_add_d_qd(qd_mul_pwr2_qd(cos(k, harf_n), 2.0), 2.0)), 0.5);    
}

QD QD::sin(int k, int n) {
    if(n <= 2)
        return 0;

    return cos(abs(k - (n >> 2)), n);
}

double QD::to_double(const QD& a) {
    return a.x[0];
}

QD& QD::operator =(const QD& r) {
    x[0] = r.x[0];
    x[1] = r.x[1];
    x[2] = r.x[2];
    x[3] = r.x[3];
    return *this;
}

QD& QD::operator =(double r) {
    x[0] = r;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    return *this;
}

QD QD::operator +(const QD& r) {
    return qd_add_qd_qd(*this, r);
}

QD QD::operator +(double r) {
    return qd_add_d_qd(*this, r);
}

QD QD_Lib::operator +(double l, const QD& r) {
    return QD::qd_add_d_qd(r, l);
}

QD QD::operator -(const QD& r) {
    return qd_sub_qd_qd(*this, r);
}

QD QD::operator -(double r) {
    return qd_add_d_qd(*this, -r);
}

QD QD_Lib::operator -(double l, const QD& r) {
    return QD::d_sub_qd_qd(l, r);
}

QD QD::operator -() {
    return QD(-x[0], -x[1], -x[2], -x[3]);
}

QD QD::operator *(const QD& r) {
    return qd_mul_qd_qd(*this, r);
}

QD QD::operator *(double r) {
    return qd_mul_d_qd(*this, r);
}

QD QD_Lib::operator *(double l, const QD& r) {
    return QD::qd_mul_d_qd(r, l);
}

QD QD::operator /(const QD& r) {
    return qd_div_qd_qd(*this, r);
}

QD QD::operator /(double r) {
    return qd_div_d_qd(*this, r);
}

QD QD_Lib::operator /(double l, const QD& r) {
    return QD::d_div_qd_qd(l, r);
}

QD QD::operator ^(int r) {
    return pow(*this, r);
}

bool QD::operator ==(const QD& r) {
    return (x[0] != r.x[0] ? 0 :
            x[1] != r.x[1] ? 0 :
            x[2] != r.x[2] ? 0 : 
            x[3] == r.x[3]);
}

bool QD::operator !=(const QD& r) {
    return (x[0] != r.x[0] ? 1 :
            x[1] != r.x[1] ? 1 :
            x[2] != r.x[2] ? 1 : 
            x[3] != r.x[3]);
}

bool QD::operator >(const QD& r) {
    return (x[0] != r.x[0] ? x[0] > r.x[0] :
            x[1] != r.x[1] ? x[1] > r.x[1] :
            x[2] != r.x[2] ? x[2] > r.x[2] : 
            x[3] > r.x[3]);
}

bool QD::operator >=(const QD& r) {
    return (x[0] != r.x[0] ? x[0] > r.x[0] :
            x[1] != r.x[1] ? x[1] > r.x[1] :
            x[2] != r.x[2] ? x[2] > r.x[2] : 
            x[3] >= r.x[3]);
}

bool QD::operator <(const QD& r) {
    return (x[0] != r.x[0] ? x[0] < r.x[0] :
            x[1] != r.x[1] ? x[1] < r.x[1] :
            x[2] != r.x[2] ? x[2] < r.x[2] : 
            x[3] < r.x[3]);
}

bool QD::operator <=(const QD& r) {
    return (x[0] != r.x[0] ? x[0] < r.x[0] :
            x[1] != r.x[1] ? x[1] < r.x[1] :
            x[2] != r.x[2] ? x[2] < r.x[2] : 
            x[3] <= r.x[3]);
}