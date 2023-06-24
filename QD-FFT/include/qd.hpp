#pragma once

namespace QD {
    using qd = double[4];

    inline void zero(qd a);
    inline void init(qd a, double x0);
    inline void init(qd a, double x0, double x1, double x2, double x3);

    void renormalize(qd a);
    void renormalize(qd a, double b);

    void add(const qd a, const qd b, qd s);
    void add(const qd a, double b, qd s);

    void sub(const qd a, const qd b, qd s);
    void sub(const qd a, double b, qd s);
    void sub(double a, const qd b, qd s);

    void mul(const qd a, const qd b, qd p);
    void mul(const qd a, double b, qd p);
    void mul(double a, double b, qd p);
    void mul_pwr2(const qd a, double b, qd p);

    void div(const qd a, const qd b, qd d);
    void div(const qd a, double b, qd d);
    void div(double a, const qd b, qd d);
    void div(double a, double b, qd d);

    void sqrt(const qd a, qd b);
}  // namespace QD