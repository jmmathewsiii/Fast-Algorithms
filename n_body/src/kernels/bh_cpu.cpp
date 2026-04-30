#include "../../include/force_cpu.h"
#include "../../include/quadtree_cpu.h"

void BH::calculate_accelerations(Stars &s, double eps, QuadTree &tree, double theta)
{
    tree.build();
    tree.compute_accelerations(eps, theta);
}
