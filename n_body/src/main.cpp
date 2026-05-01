#include "../include/initialize.h"
#include "../include/plotter.h"
#include "../include/integrator.h"
#include "../include/validate.h"
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

#if defined(BUILD_SERIAL_DIRECT)
    #define BACKEND_TAG "serial"
    #define ALGO_TAG    "direct"
    #define RUN_DIRECT  1
    #define RUN_BH      0
#elif defined(BUILD_SERIAL_BH)
    #define BACKEND_TAG "serial"
    #define ALGO_TAG    "bh"
    #define RUN_DIRECT  0
    #define RUN_BH      1
#elif defined(BUILD_CUDA_DIRECT)
    #define BACKEND_TAG "cuda"
    #define ALGO_TAG    "direct"
    #define RUN_DIRECT  1
    #define RUN_BH      0
#else
    #error "Define one of BUILD_SERIAL_DIRECT, BUILD_SERIAL_BH, BUILD_CUDA_DIRECT"
#endif

template <typename T>
static T prompt_value(const std::string& label, T default_val)
{
    std::cout << label << " [" << default_val << "]: ";
    std::string line;
    if (!std::getline(std::cin, line) || line.empty()) return default_val;
    std::istringstream iss(line);
    T val;
    if (iss >> val) return val;
    return default_val;
}

int main()
{
    const std::string backend = BACKEND_TAG;
    const std::string algo    = ALGO_TAG;
    const std::string my_tag  = algo + "_" + backend;

    std::cout << "=== N-body simulation [" << algo << " / " << backend << "] ===\n";

    int N = prompt_value<int>("Number of particles", 100);
    std::size_t n_iter = prompt_value<std::size_t>("Number of iterations", 1000);

    std::cout << "Mode (speed / plot) [plot]: ";
    std::string mode_str;
    std::getline(std::cin, mode_str);
    const bool plot = (mode_str != "speed");

    std::cout << "Running with N=" << N << ", n_iter=" << n_iter
              << ", mode=" << (plot ? "plot" : "speed") << "\n";

    double sigma_pos    = 0.3;
    double epsilon      = 0.05;
    double dt           = 0.001;
    std::uint64_t seed  = 12345;
    Random rng(seed);

    Stars stars;
    initialize(stars, N, sigma_pos, rng);

    // Snapshot dumps for Lyapunov / BH-error analysis only fire in plot mode.
    const std::size_t snapshot_stride = plot ? 10 : 0;
    const std::string snap_tag        = my_tag + "_snap";

#if RUN_DIRECT
    Direct::simulate(stars, n_iter, dt, epsilon, "direct", plot,
                     snap_tag, snapshot_stride);
#elif RUN_BH
    BH::simulate(stars, n_iter, dt, epsilon, "bh", "bh_tree", plot,
                 snap_tag, snapshot_stride);
#endif

    save_final_positions(stars, my_tag);

    // Speed-mode comparison against the direct-serial baseline (one-shot final state).
    if (!plot && my_tag != "direct_serial") {
        Stars baseline;
        if (load_final_positions(baseline, "direct_serial")) {
            compare_positions(baseline, stars, my_tag);
        } else {
            std::cout << "[validate] no direct_serial baseline found; "
                         "run 'make run-naive' (in speed mode) first.\n";
        }
    }

    // Plot-mode Lyapunov: serial-bh dumps snapshots; if serial-direct snapshots
    // also exist on disk, emit a per-step error CSV + gnuplot script.
#if RUN_BH
    if (plot && backend == "serial") {
        compare_snapshots("direct_serial_snap", "bh_serial_snap",
                          n_iter, snapshot_stride, "bh_vs_direct_serial");
    }
#endif

    return 0;
}
