#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "star.h"
#include <cstddef>
#include <string>

namespace Direct {
    void simulate(Stars& s, std::size_t n_iter, double dt, double eps,
                  const std::string& plotname, bool plot = true);
}

namespace BH {
    void simulate(Stars& s, std::size_t n_iter, double dt, double eps,
                  const std::string& plotname, const std::string& treename,
                  bool plot = true);
}

double total_energy(const Stars& s, double eps);

#endif // !INTEGRATOR_H
