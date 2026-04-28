#include <cstddef>
#include <cstdint>
#include "star.h"

constexpr idx NO_CHILD = UINT32_MAX;

struct Node
{
    /* Always stored */
    double cx, cy;
    double mass;
    double size;

    /* Index offset into the tree */
    idx first_child;

    /* Index offset into Stars struct-of-arrays */
    idx star_idx;

};

