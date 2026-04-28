#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "random.h"
#include "star.h"

void initialize(Stars& stars, int N, double sigma_pos, Random& rng);

void virial_rescale(Stars& stars, double eps);

#endif // !INITIALIZE_H
