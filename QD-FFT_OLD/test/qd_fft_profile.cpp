#include <gperftools/profiler.h>

#include "fft.hpp"
#include "qd.hpp"

int main() {
    const uint16_t N = (1 << 15);

    QD::qd cos_table[N / 4 + 1];
    QD::qd x_[N];
    QD::qd ix_[N];
    QD::qd y_[N];
    QD::qd iy_[N];
    double* x[N];
    double* ix[N];
    double* y[N];
    double* iy[N];
    FFT::make_cos_table(N, cos_table);

    for (uint16_t j = 0; j < N; j++) {
        QD::rand(x_[j]);
        QD::rand(ix_[j]);
        x[j]  = x_[j];
        ix[j] = ix_[j];
        y[j]  = y_[j];
        iy[j] = iy_[j];
    }
    ProfilerStart("qd_fft_profile.prof");
    FFT::fft(N, N, x, ix, y, iy, cos_table);
    ProfilerStop();
}