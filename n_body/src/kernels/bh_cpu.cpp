#include "../../include/force_cpu.h"
#include "../../include/quadtree_cpu.h"
#include <chrono>

namespace {
    double g_build_seconds = 0.0;
    double g_force_seconds = 0.0;
}

void BH::reset_timers()
{
    g_build_seconds = 0.0;
    g_force_seconds = 0.0;
}

double BH::tree_build_seconds()   { return g_build_seconds; }
double BH::force_compute_seconds() { return g_force_seconds; }

void BH::calculate_accelerations(Stars &s, double eps, QuadTree &tree, double theta)
{
    auto t0 = std::chrono::steady_clock::now();
    tree.build();
    auto t1 = std::chrono::steady_clock::now();
    tree.compute_accelerations(eps, theta);
    auto t2 = std::chrono::steady_clock::now();

    g_build_seconds += std::chrono::duration<double>(t1 - t0).count();
    g_force_seconds += std::chrono::duration<double>(t2 - t1).count();
}
