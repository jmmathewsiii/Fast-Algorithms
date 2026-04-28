#include <cstdint>
#include <iostream>
#include <cmath>
#include <chrono>
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

    // hodlr.printMetaData();
    // hodlr.printTreeData();

    VD x(N);
    for (std::size_t i = 0; i < N; ++i) {
        x[i] = rng.uniform();
    }
    auto hodlr_start = std::chrono::high_resolution_clock::now();
    VD b = hodlr.MatVec(x);
    auto hodlr_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> hodlr_dur = hodlr_end - hodlr_start;

    VD A = hodlr.createFullMatrix();
    VD bp(N);
    
    auto normal_start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            bp[i] += A[(i * N) + j] * x[j];
    auto normal_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> normal_dur = normal_end - normal_start;

    VD err;
    err.reserve(N);
    for (std::size_t i = 0; i < N; ++i)
        err.push_back(std::abs(b[i] - bp[i]));

    double sum = 0;
    for (std::size_t i = 0; i < N; ++i)
        sum += err[i];

    std::cout << "Sum of errs: " << sum << "\n";
    std::cout << "HODLR MatVec Duration: " << hodlr_dur.count() << " ms.\n";
    std::cout << "Normal MatVec Duration: " << normal_dur.count() << " ms.\n";

    return 0;
}

