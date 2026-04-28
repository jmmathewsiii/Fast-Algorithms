#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "../include/star.h"

void calculate_accelerations(Stars&, double);

void kdk_verlet_step(Stars&, double, double);

double total_energy(const Stars&, double);

void simulate(Stars&, std::size_t, double, double);


#endif // !INTEGRATOR_H
