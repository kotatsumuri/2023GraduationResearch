#pragma once
#include <random>

#include "qd.hpp"

void rand(qd a) {
    static std::random_device rnd;
    static std::mt19937 engine{rnd()};
    static std::uniform_real_distribution<> dist{0, 1};
    static const double m_const = 4.6566128730773926e-10;
    double m                    = 1;
    zero(a);
    double d;

    for (int i = 0; i < 7; i++, m *= m_const) {
        d = dist(engine) * m;
        add(a, d, a);
    }
}