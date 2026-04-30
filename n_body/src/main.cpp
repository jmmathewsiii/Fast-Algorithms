#include "../include/initialize.h"
#include "../include/plotter.h"
#include "../include/integrator.h"
#include <cstdint>
#include <iostream>

int main()
{
    Stars stars;

    int N = 10000;
    std::size_t n_iter = 4000;
    double sigma_pos = 0.3;
    double epsilon = 0.05;
    double dt = 0.001;
    std::uint64_t seed = 12345;
    Random rng(seed);

    initialize(stars, N, sigma_pos, rng);
    // stars.m.push_back(10.);
    // stars.x.push_back(3.);
    // stars.y.push_back(3.);
    // stars.vx.push_back(0.);
    // stars.vy.push_back(0.);
    // stars.ax.push_back(0.);
    // stars.ay.push_back(0.);

    Plotter::plot_positions(stars, "initial");

    // Snapshot initial conditions so both runs see the same start.
    Stars stars_initial = stars;

    // std::cout << "=== Direct O(N^2) run ===\n";
    // Direct::simulate(stars, n_iter, dt, epsilon, "direct");

    // Restore initial conditions for the BH run.
    stars = stars_initial;

    std::cout << "=== Barnes-Hut run ===\n";
    BH::simulate(stars, n_iter, dt, epsilon, "bh", "bh_tree");

    return 0;
}
