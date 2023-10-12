#pragma once

namespace util::calc {
void split(double a, double* a_hi, double* a_lo);
void quick_two_sum(double a, double b, double* s, double* e);
void two_sum(double a, double b, double* s, double* e);
void three_sum(double a, double b, double c, double* s, double* e0, double* e1);
void three_sum(double a, double b, double c, double* s, double* e0);
void two_prod(double a, double b, double* p, double* e);
void six_three_sum(double a, double b, double c, double d, double e, double f,
                   double* s, double* e0, double* e1);
void dd_add_dd_dd(double a0, double a1, double b0, double b1, double* s0,
                  double* s1);
void nine_two_sum(double a, double b, double c, double d, double e, double f,
                  double g, double h, double i, double* s, double* e0);

int ufp(double a);
int ulp(double a);
int sign(double a);
}  // namespace util::calc