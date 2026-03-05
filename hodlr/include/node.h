#ifndef NODE_H
#define NODE_H

#include <cstddef>

struct Node {

    /* Always stored */
    std::size_t level; // Level 0 is full structure; L is smallest (leaves)
    std::size_t start; // Where this node would start in the full N x N matrix
    std::size_t end; // Where this node would end in the full N x N matrix (noninclusive)
    std::size_t length = start - end;
    std::size_t uv_length = length / static_cast<std::size_t>(2);
    bool isLeaf; // Is it a leaf node (dense n x n matrix)

    /* If NOT a leaf */
    std::size_t left_child; // Index of A11 in the node array
    std::size_t right_child; // Index of A22 in the node array
    std::size_t U1; // Index of start of U12 in U array
    std::size_t U2; // Index of start of U21 in U array
    std::size_t VT1; // Index of start of VT12 in V transpose array
    std::size_t VT2; // Index of start of VT21 in V transpose array

    /* If leaf */
    std::size_t A; // Index of start of A in the leaf array

};

#endif // !NODE_H
