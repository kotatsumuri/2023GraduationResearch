#include <cstdio>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <qd/qd_real.h>
#include "util.hpp"
#include "qd.hpp"

using namespace QD_Lib;

namespace mp = boost::multiprecision;
typedef mp::number<mp::cpp_bin_float<1000>> Real;
Real mp_pi = mp::acos(static_cast<Real>(-1.0));

Real reccos(unsigned long long int k, unsigned long long n) {
    if(k == 0)
        return static_cast<Real>(1.0);
    if(n == 1)
        return static_cast<Real>(1.0);
    if(n == 2)
        return static_cast<Real>(-1.0);

    unsigned long long int harf_n = n >> 1;

    if(k > harf_n)
        k = n - k;
    
    if(k > (harf_n >> 1))
        return - mp::sqrt(2.0 + 2.0 * reccos(harf_n - k, harf_n)) / 2.0;
    else 
        return mp::sqrt(2.0 + 2.0 * reccos(k, harf_n)) / 2.0;
}

int main() {
    unsigned long long int n = 1;
    for(int i = 0;i < 64;i++) {
        Real a = reccos(1, n);
        QD b = QD::cos(1, n);
        Real e = mp::abs(a - b.x[0] - b.x[1] - b.x[2] - b.x[3]);
        std::cout << n << "," << e.str(0, std::ios_base::scientific) << std::endl;
        n *= 2;
    }

    return 0;
}
