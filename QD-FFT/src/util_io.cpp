#include "util_io.hpp"

#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "util_calc.hpp"

namespace util::io {
std::string to_bin_string(double a) {
    std::string ret = "";

    if (a == 0) {
        ret = "zero";
        return ret;
    }

    if (std::isinf(a)) {
        ret = "inf";
        return ret;
    }

    if (std::isnan(a)) {
        ret = "nan";
        return ret;
    }

    if (util::calc::sign(a) > 0)
        ret += "+2^";
    else
        ret += "-2^";

    if (util::calc::ufp(a) >= 0)
        ret += "+";
    else
        ret += "-";

    if (util::calc::ufp(a) == -1023) {
        ret += "1022";
        ret += "*";
    } else {
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << std::abs(util::calc::ufp(a));
        ret += ss.str() + "*";
    }
    if (util::calc::ufp(a) == -1023)
        ret += "0.";
    else
        ret += "1.";

    unsigned long long row_bit = *(unsigned long long *)&a;
    std::bitset<64> bin(row_bit);

    ret += bin.to_string().substr(12);

    return ret;
}
}  // namespace util::io