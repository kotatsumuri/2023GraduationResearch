#pragma once
#include <string>

namespace util {
    std::string to_reg_str(double a);
    void split(double a, double* a_hi, double* a_lo);
    void quick_two_sum(double a, double b, double* s, double* e);
    void two_sum(double a, double b, double* s, double* e);
    void three_sum(double a, double b, double c, double* s, double* e0, double* e1);
    void three_sum(double a, double b, double c, double* s, double* e0);
    void two_prod(double a, double b, double* p, double* e);
}