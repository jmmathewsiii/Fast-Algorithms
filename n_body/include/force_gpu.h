#ifndef FORCE_GPU_H
#define FORCE_GPU_H

#include "star.h"

class QuadTree;

namespace Direct {
    void calculate_accelerations(DeviceStars s, double eps);
}

namespace BH {
    void calculate_accelerations(DeviceStars s, double eps, QuadTree& tree, double theta);
}

#endif // !FORCE_GPU_H
