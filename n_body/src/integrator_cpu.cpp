#include "../include/integrator.h"
#include "../include/plotter.h"
#include "../include/force_cpu.h"
#include "../include/quadtree_cpu.h"
#include "../include/validate.h"
#include <chrono>
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
                      const std::string& plotname, bool plot,
                      const std::string& snapshot_tag,
                      std::size_t snapshot_stride)
{
    Direct::reset_timers();
    auto t_start = std::chrono::steady_clock::now();

    Direct::calculate_accelerations(s, eps);

    Plotter::Animator anim;
    if (plot) {
        anim = Plotter::animator_begin(plotname);
        Plotter::animator_add_frame(anim, s);
    }

    const bool snap = snapshot_stride > 0 && !snapshot_tag.empty();
    if (snap) save_snapshot(s, snapshot_tag, 0);

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kick_drift(s, dt);
        Direct::calculate_accelerations(s, eps);
        kick(s, dt);

        if (plot && i % 10 == 9)
        {
            Plotter::animator_add_frame(anim, s);
        }
        if (!plot && i % 200 == 199)
        {
            double energy = total_energy(s, eps);
            std::cout << "[direct] Total Energy at step " << i + 1 << ": " << energy << "\n";
        }
        if (snap && (i + 1) % snapshot_stride == 0)
        {
            save_snapshot(s, snapshot_tag, i + 1);
        }
    }

    if (plot) {
        Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
    }

    if (!plot) {
        double total = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();
        double force = Direct::force_compute_seconds();
        std::cout << "[direct] Total run time:  " << total << " s\n"
                  << "[direct] Force calc time: " << force << " s\n";
    }
}

void BH::simulate(Stars &s, std::size_t n_iter, double dt, double eps,
                  const std::string& plotname, const std::string& treename,
                  bool plot,
                  const std::string& snapshot_tag,
                  std::size_t snapshot_stride)
{
    constexpr double theta = 0.5;
    QuadTree tree(&s);

    BH::reset_timers();
    auto t_start = std::chrono::steady_clock::now();

    BH::calculate_accelerations(s, eps, tree, theta);

    Plotter::Animator anim;
    Plotter::TreeAnimator tree_anim;
    if (plot) {
        anim = Plotter::animator_begin(plotname);
        tree_anim = Plotter::tree_animator_begin(treename);
        Plotter::animator_add_frame(anim, s);
        Plotter::tree_animator_add_frame(tree_anim, tree, s);
    }

    const bool snap = snapshot_stride > 0 && !snapshot_tag.empty();
    if (snap) save_snapshot(s, snapshot_tag, 0);

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kick_drift(s, dt);
        BH::calculate_accelerations(s, eps, tree, theta);
        kick(s, dt);

        if (plot && i % 10 == 9)
        {
            Plotter::animator_add_frame(anim, s);
            Plotter::tree_animator_add_frame(tree_anim, tree, s);
        }
        if (!plot && i % 200 == 199)
        {
            double energy = total_energy(s, eps);
            std::cout << "[bh]     Total Energy at step " << i + 1 << ": " << energy << "\n";
        }
        if (snap && (i + 1) % snapshot_stride == 0)
        {
            save_snapshot(s, snapshot_tag, i + 1);
        }
    }

    if (plot) {
        Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
        Plotter::tree_animator_end(tree_anim, /*half_range=*/4.0, /*pause_sec=*/0.005);
    }

    if (!plot) {
        double total = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();
        double build = BH::tree_build_seconds();
        double force = BH::force_compute_seconds();
        std::cout << "[bh]     Total run time:   " << total << " s\n"
                  << "[bh]     Tree build time:  " << build << " s\n"
                  << "[bh]     Force calc time:  " << force << " s\n";
    }
}
