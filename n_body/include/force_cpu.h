#ifndef FORCE_CPU_H
#define FORCE_CPU_H

#include "star.h"

class QuadTree;

namespace Direct {
    void calculate_accelerations(Stars& s, double eps);
}

namespace BH {
    void calculate_accelerations(Stars& s, double eps, QuadTree& tree, double theta);

    // Timing accumulators, in seconds. Reset before a run; query after.
    void reset_timers();
    double tree_build_seconds();
    double force_compute_seconds();
}

#endif // !FORCE_CPU_H
