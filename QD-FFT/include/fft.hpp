#pragma once
#include "qd.hpp"

namespace FFT {
void make_cos_table(long long int N, QD::qd cos_table[]);
void decimation_in_frequency(long long int initN, long long int n, QD::qd x[],
                             QD::qd ix[], QD::qd y[], QD::qd iy[],
                             QD::qd cos_table[]);
void inv_decimation_in_frequency(long long int initN, long long int n,
                                 QD::qd x[], QD::qd ix[], QD::qd y[],
                                 QD::qd iy[], QD::qd cos_table[]);
void DFT(int N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
         QD::qd cos_table[]);
void inv_DFT(int N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
             QD::qd cos_table[]);
}  // namespace FFT