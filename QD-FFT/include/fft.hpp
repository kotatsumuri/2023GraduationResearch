#pragma once
#include "qd.hpp"

namespace FFT {
void make_cos_table(uint16_t N, QD::qd cos_table[]);
void make_cos_table(uint16_t N, double cos_table[]);
void decimation_in_frequency(uint16_t initN, uint16_t n, QD::qd x[],
                             QD::qd ix[], QD::qd y[], QD::qd iy[],
                             QD::qd cos_table[]);
void decimation_in_frequency(uint16_t initN, uint16_t n, double x[],
                             double ix[], double y[], double iy[],
                             double cos_table[]);
void inv_decimation_in_frequency(uint16_t initN, uint16_t n, QD::qd x[],
                                 QD::qd ix[], QD::qd y[], QD::qd iy[],
                                 QD::qd cos_table[]);
void inv_decimation_in_frequency(uint16_t initN, uint16_t n, double x[],
                                 double ix[], double y[], double iy[],
                                 double cos_table[]);
void DFT(uint16_t N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
         QD::qd cos_table[]);
void inv_DFT(uint16_t N, QD::qd x[], QD::qd ix[], QD::qd y[], QD::qd iy[],
             QD::qd cos_table[]);
}  // namespace FFT