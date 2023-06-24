#pragma once

namespace QD {
    using qd = double[4];

    void renormalize(qd a);
    void renormalize(qd a, double b);

    void add(const qd a, const qd b, qd s);
    void add(const qd a, double b, qd s);
}  // namespace QD