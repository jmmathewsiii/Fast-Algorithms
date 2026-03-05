#include <cstdint>
#include <iostream>
#include <cmath>
#include "../include/hodlr.h"

void printMetaData(HODLR_Matrix&);

int main(int argc, char** argv)
{

    std::size_t n = static_cast<std::size_t>(std::stoull(argv[1]));        // Size of leaf nodes: n x n
    std::size_t k = static_cast<std::size_t>(std::stoull(argv[2]));        // Rank of low-rank blocks
    std::size_t L = static_cast<std::size_t>(std::stoull(argv[3]));        // Level of leaf nodes
    std::uint64_t seed = static_cast<std::uint64_t>(std::stoull(argv[4])); // Random seed

    std::size_t N = (static_cast<std::size_t>(1) << L) * n;

    Random rng(seed);
    HODLR_Matrix hodlr(n, k, L);
    hodlr.fillWithRandomData(rng);
    hodlr.initializeTree();

    hodlr.printMetaData();
    hodlr.printTreeData();

    VD x(N);
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = rng.uniform();
    }
    // VD b = hodlr.MatVec(x);

     return 0;
}

