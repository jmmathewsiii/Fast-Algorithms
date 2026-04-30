#ifndef PLOTTER_H
#define PLOTTER_H

#include <string>
#include "star.h"

class QuadTree;

namespace Plotter {
    void plot_positions(const Stars& stars, const std::string& name);

    struct Animator {
        std::string name;
        std::size_t n_frames = 0;
    };

    Animator animator_begin(const std::string& name);
    void     animator_add_frame(Animator& a, const Stars& s);
    void     animator_end(Animator& a, double half_range, double pause_sec);

    struct TreeAnimator {
        std::string name;
        std::size_t n_frames = 0;
    };

    TreeAnimator tree_animator_begin(const std::string& name);
    void         tree_animator_add_frame(TreeAnimator& a, const QuadTree& tree, const Stars& s);
    void         tree_animator_end(TreeAnimator& a, double half_range, double pause_sec);
}

#endif // !PLOTTER_H
