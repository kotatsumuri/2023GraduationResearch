#pragma once

namespace QD {
    using qd = double[4];

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
}  // namespace QD