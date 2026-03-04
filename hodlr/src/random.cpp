#include "../include/random.h"

Random::Random(std::uint64_t seed)
    : generator(seed),
      u_dist(0.0, 1.0),
      n_dist(0.0, 1.0)
{}

double Random::uniform()
{
    return u_dist(generator);
}

double Random::normal()
{
    return n_dist(generator);
}
