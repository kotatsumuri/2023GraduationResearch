#pragma once

using qd = double[4];

struct qd_complex {
    qd re;
    qd im;
};

#include "common.hpp"
#include "qd_base.hpp"
#include "qd_init.hpp"
#include "qd_io.hpp"
#include "qd_random.hpp"
#include "qd_sqrt.hpp"
#include "qd_util.hpp"
#include "qd_vector.hpp"