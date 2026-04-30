#ifndef FORCE_CPU_H
#define FORCE_CPU_H

#include "star.h"

class QuadTree;

namespace Direct {
    void calculate_accelerations(Stars& s, double eps);
}

namespace BH {
    void calculate_accelerations(Stars& s, double eps, QuadTree& tree, double theta);
}

#endif // !FORCE_CPU_H
