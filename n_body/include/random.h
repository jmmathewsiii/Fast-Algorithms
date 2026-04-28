#ifndef RANDOM_H
#define RANDOM_H

#include <cstdint>
#include <random>

class Random {
public:
    explicit Random(std::uint64_t seed);

    double uniform();
    double normal();

private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> u_dist;
    std::normal_distribution<double> n_dist;
};

#endif // !RANDOM_H
