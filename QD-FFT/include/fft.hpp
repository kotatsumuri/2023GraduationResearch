#pragma once
#include "qd.hpp"

namespace FFT {
void make_cos_table(int N, QD::qd cos_table[]);
void decimation_in_frequency(int initN, int n, QD::qd x[], QD::qd ix[],
                             QD::qd y[], QD::qd iy[], QD::qd cos_table[]);
void inv_decimation_in_frequency(int initN, int n, QD::qd x[], QD::qd ix[],
                                 QD::qd y[], QD::qd iy[], QD::qd cos_table[]);
}  // namespace FFT