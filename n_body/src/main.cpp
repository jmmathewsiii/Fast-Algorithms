#include "../include/initialize.h"
#include "../include/plotter.h"
#include "../include/integrator.h"
#include "../include/validate.h"
#include <cstdint>
#include <iostream>
#include <string>

#ifndef BACKEND_TAG
#define BACKEND_TAG "unknown"
#endif

int main()
{
    const std::string backend = BACKEND_TAG;

    // Set to false for clean timing runs (no per-frame file I/O).
    const bool plot = false;

    Stars stars;

    int N = 10000;
    std::size_t n_iter = 4000;
    double sigma_pos = 0.3;
    double epsilon = 0.05;
    double dt = 0.001;
    std::uint64_t seed = 12345;
    Random rng(seed);

    initialize(stars, N, sigma_pos, rng);

    if (plot) Plotter::plot_positions(stars, "initial");

    // Snapshot initial conditions so all runs see the same start.
    Stars stars_initial = stars;

    std::cout << "=== Direct O(N^2) run [" << backend << "] ===\n";
    Direct::simulate(stars, n_iter, dt, epsilon, "direct", plot);
    save_final_positions(stars, "direct_" + backend);

    // If this is the serial run, this IS the baseline; otherwise compare to it.
    if (backend != "serial") {
        Stars baseline;
        if (load_final_positions(baseline, "direct_serial")) {
            compare_positions(baseline, stars, "direct-" + backend);
        } else {
            std::cout << "[validate] no direct_serial baseline found; "
                         "run the serial binary first.\n";
        }
    }

    // Restore initial conditions for the BH run.
    stars = stars_initial;

    std::cout << "=== Barnes-Hut run [" << backend << "] ===\n";
    BH::simulate(stars, n_iter, dt, epsilon, "bh", "bh_tree", plot);
    save_final_positions(stars, "bh_" + backend);

    {
        Stars baseline;
        if (load_final_positions(baseline, "direct_serial")) {
            compare_positions(baseline, stars, "bh-" + backend);
        } else {
            std::cout << "[validate] no direct_serial baseline found; "
                         "run the serial binary first.\n";
        }
    }

    return 0;
}
