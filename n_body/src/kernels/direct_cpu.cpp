#include "../../include/force_cpu.h"
#include "../../include/star.h"
#include <chrono>
#include <cmath>

namespace {
    double g_force_seconds = 0.0;
}

void Direct::reset_timers() { g_force_seconds = 0.0; }
double Direct::force_compute_seconds() { return g_force_seconds; }

void Direct::calculate_accelerations(Stars &s, double eps)
{
    auto t0 = std::chrono::steady_clock::now();

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

    g_force_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}
