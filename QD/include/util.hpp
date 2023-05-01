#pragma once
#include <string>

namespace util {
    std::string to_reg_str(double a);
    void split(double a, double* a_hi, double* a_lo);
    double quick_two_sum(double a, double b, double* e);
    double two_sum(double a, double b, double* e);
    double three_sum(double a, double b, double c, double* e0, double* e1);
    double three_sum(double a, double b, double c, double* e0);
    double two_prod(double a, double b, double* e);
    double six_three_sum(double a, double b, double c, double d, double e, double f, double* e0, double* e1);
    void dd_add_dd_dd(double a0, double a1, double b0, double b1, double* s0, double* s1);
    double nine_two_sum(double a, double b, double c, double d, double e, double f, double g, double h, double i, double* e0);
}