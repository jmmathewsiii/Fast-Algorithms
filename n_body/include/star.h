#ifndef STAR_H
#define STAR_H

#include <cstdint>
#include <vector>

using VD = std::vector<double>;
using idx = std::uint32_t;

constexpr double G = 1.;

struct Stars {
    VD m;
    VD x, y;
    VD vx, vy;
    VD ax, ay;

    idx size() const { return x.size(); }

    void resize(idx N)
    {
        m.resize(N);
        x.resize(N); y.resize(N);
        vx.resize(N); vy.resize(N);
        ax.resize(N); ay.resize(N);
    }
};


#endif // !STAR_H
