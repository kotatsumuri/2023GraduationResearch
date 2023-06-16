#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <bitset>
#include "QD.hpp"

using namespace QD_Lib;

namespace mp = boost::multiprecision;
const int BITS = 1285;
typedef mp::number<mp::cpp_bin_float<BITS, mp::backends::digit_base_2>> Real;
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

int ulp(double a) {
    unsigned long long row_bit = *(unsigned long long *)&a;
    std::bitset<64> bin(row_bit);
    return ((int)(row_bit >> 52) & 2047) - 1023 - 52;
}

int main(void) {
    std::cout << "static const QD cos_table[] = {" << std::endl;

    unsigned long long n = 2048;
    for(unsigned long long k = 0;k <= n / 4;k++) {
        QD a(0.0);
        Real b = reccos(k, n);
        int ulp_b = static_cast<int>(mp::floor(mp::log2(b))) + 1 - (BITS - 1);
        Real c = mp::ldexp(b, -ulp_b);
        if(b != 0) {
            for(int i = 0;i < BITS;i++) {
                a = QD::qd_add_d_qd(a, static_cast<double>(mp::ldexp(c - mp::ldexp(mp::floor(mp::ldexp(c, -1)), 1),ulp_b)));
                c = mp::floor(mp::ldexp(c, -1));
                ulp_b++;
            }
        }
        // C++コードを生成 
        std::cout << "    QD(";
        printf("%.16e, %.16e,\n", a.x[0], a.x[1]);
        printf("       %.16e, %.16e)", a.x[2], a.x[3]);
        if(k != n / 4)
            std::cout << ",";
        std::cout << std::endl;
        
        // 誤差を計算
        // Real e = mp::abs(b - a.x[0] - a.x[1] - a.x[2] - a.x[3]);
        // Real error_bits = 0;
        // if(e != 0)
        //     error_bits =  mp::log2(e) - static_cast<Real>(ulp(a.x[3]));
        // std::cout << k << "," << error_bits.str(0, std::ios_base::scientific) << std::endl;
    }

    std::cout << "};" << std::endl;
}