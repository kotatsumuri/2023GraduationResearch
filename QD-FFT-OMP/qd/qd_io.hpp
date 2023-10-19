#pragma once
#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "qd.hpp"
#include "qd_util.hpp"

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

    if (a >= 0)
        ret += "+2^";
    else
        ret += "-2^";

    if (log_ufp(a) >= 0)
        ret += "+";
    else
        ret += "-";

    if (log_ufp(a) == -1023) {
        ret += "1023";
        ret += "*";
    } else {
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << std::abs(log_ufp(a));
        ret += ss.str() + "*";
    }
    if (log_ufp(a) == -1023)
        ret += "0.";
    else
        ret += "1.";

    unsigned long long row_bit = *(unsigned long long *)&a;
    std::bitset<64> bin(row_bit);

    ret += bin.to_string().substr(12);

    return ret;
}

std::string to_bin_string(const qd a) {
    return to_bin_string(a[0]) + "|" + to_bin_string(a[1]) + "|" +
           to_bin_string(a[2]) + "|" + to_bin_string(a[3]);
}