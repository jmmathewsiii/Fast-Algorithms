#ifndef PLOTTER_H
#define PLOTTER_H

#include <fstream>
#include <string>
#include "star.h"

namespace Plotter {
    void plot_positions(const Stars& stars, const std::string& name);

    struct Animator {
        std::ofstream data;
        std::string   name;
        std::size_t   n_frames = 0;
    };

    Animator animator_begin(const std::string& name);
    void     animator_add_frame(Animator& a, const Stars& s);
    void     animator_end(Animator& a, double half_range, double pause_sec);
}

#endif // !PLOTTER_H
