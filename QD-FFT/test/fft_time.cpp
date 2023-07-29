#include <chrono>
#include <iostream>
#include <random>

#include "fft.hpp"
#include "qd.hpp"

int main() {
    const uint16_t M = 15;
    const uint16_t K = 10;
    uint16_t N       = 2;
    for (uint16_t i = 0; i < M; i++, N *= 2) {
        double_t microsec = 0;
        for (uint16_t k = 0; k < K; k++) {
            QD::qd cos_table[N / 4 + 1];
            QD::qd x[N];
            QD::qd ix[N];
            QD::qd y[N];
            QD::qd iy[N];

            FFT::make_cos_table(N, cos_table);

            for (uint16_t j = 0; j < N; j++) {
                QD::rand(x[j]);
                QD::rand(ix[j]);
            }

            auto start = std::chrono::system_clock::now();
            FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);
            auto end = std::chrono::system_clock::now();

            auto dur = end - start;
            microsec +=
                std::chrono::duration_cast<std::chrono::microseconds>(dur)
                    .count();
        }
        std::cout << N << ", " << std::scientific << microsec / K << std::endl;
    }

    std::random_device rnd;
    std::mt19937 engine{rnd()};
    std::uniform_real_distribution<> dist{0, 1};
    N = 2;
    for (uint16_t i = 0; i < M; i++, N *= 2) {
        double_t microsec = 0;
        for (uint16_t k = 0; k < K; k++) {
            double cos_table[N / 4 + 1];
            double x[N];
            double ix[N];
            double y[N];
            double iy[N];

            FFT::make_cos_table(N, cos_table);

            for (uint16_t j = 0; j < N; j++) {
                x[j]  = dist(engine);
                ix[j] = dist(engine);
            }

            auto start = std::chrono::system_clock::now();
            FFT::decimation_in_frequency(N, N, x, ix, y, iy, cos_table);
            auto end = std::chrono::system_clock::now();

            auto dur = end - start;
            microsec +=
                std::chrono::duration_cast<std::chrono::microseconds>(dur)
                    .count();
        }
        std::cout << N << ", " << std::scientific << microsec / K << std::endl;
    }
}