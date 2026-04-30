#include "../include/integrator.h"
#include <cmath>

double total_energy(const Stars &s, double eps)
{
    const std::size_t N = s.size();
    const double eps_sq = eps * eps;

    double K = 0.;
    for (std::size_t i = 0; i < N; ++i)
    {
        K += 0.5 * s.m[i] * (s.vx[i] * s.vx[i] + s.vy[i] * s.vy[i]);
    }

    double U = 0.;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
        {
            double dx = s.x[i] - s.x[j];
            double dy = s.y[i] - s.y[j];
            double r  = std::sqrt(dx * dx + dy * dy + eps_sq);
            U -= G * s.m[i] * s.m[j] / r;
        }

    return K + U;
}
