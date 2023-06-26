#include <string>

#include "qd.hpp"
#include "util_io.hpp"

namespace QD {
std::string to_bin_string(qd a) {
    return util::io::to_bin_string(a[0]) + "|" + util::io::to_bin_string(a[1]) +
           "|" + util::io::to_bin_string(a[2]) + "|" +
           util::io::to_bin_string(a[3]);
}
}  // namespace QD