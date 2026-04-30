#include "../include/integrator.h"
#include "../include/plotter.h"
#include "../include/force_cpu.h"
#include "../include/quadtree_cpu.h"
#include <iostream>
#include <cmath>

namespace {

void kick_drift(Stars &s, double dt)
{
    std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i)
    {
        s.vx[i] += 0.5 * s.ax[i] * dt;
        s.vy[i] += 0.5 * s.ay[i] * dt;
        s.x[i]  += s.vx[i] * dt;
        s.y[i]  += s.vy[i] * dt;
    }
}

void kick(Stars &s, double dt)
{
    std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i)
    {
        s.vx[i] += 0.5 * s.ax[i] * dt;
        s.vy[i] += 0.5 * s.ay[i] * dt;
    }
}

}  // namespace

void Direct::simulate(Stars &s, std::size_t n_iter, double dt, double eps,
                      const std::string& plotname)
{
    Direct::calculate_accelerations(s, eps);

    Plotter::Animator anim = Plotter::animator_begin(plotname);
    Plotter::animator_add_frame(anim, s);

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kick_drift(s, dt);
        Direct::calculate_accelerations(s, eps);
        kick(s, dt);

        if (i % 10 == 9)
        {
            Plotter::animator_add_frame(anim, s);
        }
        if (i % 200 == 199)
        {
            double energy = total_energy(s, eps);
            std::cout << "[direct] Total Energy at step " << i + 1 << ": " << energy << "\n";
        }
    }

    Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
}

void BH::simulate(Stars &s, std::size_t n_iter, double dt, double eps,
                  const std::string& plotname, const std::string& treename)
{
    constexpr double theta = 0.5;
    QuadTree tree(&s);

    BH::calculate_accelerations(s, eps, tree, theta);

    Plotter::Animator anim = Plotter::animator_begin(plotname);
    Plotter::TreeAnimator tree_anim = Plotter::tree_animator_begin(treename);
    Plotter::animator_add_frame(anim, s);
    Plotter::tree_animator_add_frame(tree_anim, tree, s);

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kick_drift(s, dt);
        BH::calculate_accelerations(s, eps, tree, theta);
        kick(s, dt);

        if (i % 10 == 9)
        {
            Plotter::animator_add_frame(anim, s);
            Plotter::tree_animator_add_frame(tree_anim, tree, s);
        }
        if (i % 200 == 199)
        {
            double energy = total_energy(s, eps);
            std::cout << "[bh]     Total Energy at step " << i + 1 << ": " << energy << "\n";
        }
    }

    Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
    Plotter::tree_animator_end(tree_anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
}
