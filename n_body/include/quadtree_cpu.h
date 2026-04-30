#ifndef QUADTREE_H
#define QUADTREE_H
#include "node.h"
#include "star.h"
#include <cstdint>
#include <iosfwd>

using Tree = std::vector<Node>;
using idx_list = std::vector<idx>;


class QuadTree
{
    private:
        Stars *s;
        Tree tree;

        double root_x_min, root_y_min, root_box_size;

        const std::uint16_t MAX_DEPTH = 50;

        void build_recursive(idx curr_idx, idx_list &star_idx_lst, double x_min, double y_min, double box_size, std::uint16_t curr_depth);
        void compute_acc_for(idx node_idx, idx star_idx, double &ax, double &ay, double eps, double theta) const;
        void write_divisions_recursive(idx node_idx, double x_min, double y_min, double box_size, std::ostream& os) const;

    public:
        QuadTree(Stars* s_);
        void build();
        void compute_accelerations(double eps, double theta) const;
        void write_divisions(std::ostream& os) const;
        const Tree& getTree() const;
};

#endif // !QUADTREE_H
