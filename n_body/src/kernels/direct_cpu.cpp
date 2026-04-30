#include "../../include/force_cpu.h"
#include "../../include/star.h"
#include <cmath>

void Direct::calculate_accelerations(Stars &s, double eps)
{
    std::size_t N = s.size();
    double eps_sq = eps * eps;

    for (std::size_t i = 0; i < N; ++i)
    {
        s.ax[i] = 0.;
        s.ay[i] = 0.;
    }

    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
        {
            double dx = s.x[i] - s.x[j];
            double dy = s.y[i] - s.y[j];
            double r_sq = dx * dx + dy * dy + eps_sq;
            double r3_inv = 1. / (r_sq * std::sqrt(r_sq));
            double ax_new = G * dx * r3_inv;
            double ay_new = G * dy * r3_inv;
            s.ax[i] -= s.m[j] * ax_new;    s.ay[i] -= s.m[j] * ay_new;
            s.ax[j] += s.m[i] * ax_new;    s.ay[j] += s.m[i] * ay_new;
        }
}
