#include "../include/integrator.h"
#include "../include/plotter.h"
#include "../include/force_cpu.h"
#include <iostream>
#include <cmath>

void kdk_verlet_step(Stars &s, double dt, double eps)
{
    std::size_t N = s.size();

    for (std::size_t i = 0; i < N; ++i)
    {
        // Kick
        s.vx[i] += 0.5 * s.ax[i] * dt;
        s.vy[i] += 0.5 * s.ay[i] * dt;

        // Drift
        s.x[i] += s.vx[i] * dt;
        s.y[i] += s.vy[i] * dt;
    }

    calculate_accelerations(s, eps);

    // Kick (again)
    for (std::size_t i = 0; i < N; ++i)
    {
        s.vx[i] += 0.5 * s.ax[i] * dt;
        s.vy[i] += 0.5 * s.ay[i] * dt;
    }
}

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

void simulate(Stars &s, std::size_t n_iter, double dt, double eps)
{
    calculate_accelerations(s, eps);

    Plotter::Animator anim = Plotter::animator_begin("nbody");
    Plotter::animator_add_frame(anim, s);

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kdk_verlet_step(s, dt, eps);
        if (i % 10 == 9)
        {
            Plotter::animator_add_frame(anim, s);
        }
        if (i % 200 == 199)
        {
            double energy = total_energy(s, eps);
            std::cout << "Total Energy at step " << i + 1 << ": " << energy << "\n";
        }
    }

    Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
}
