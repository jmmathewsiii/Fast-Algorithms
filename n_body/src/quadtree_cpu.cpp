#include "../include/quadtree_cpu.h"
#include <algorithm>
#include <cmath>
#include <cstdint>

QuadTree::QuadTree(Stars *s_) { s = s_; }

const Tree& QuadTree::getTree() const { return tree; }

void QuadTree::build()
{
    idx N = s->size();
    idx_list star_idx_lst;

    double x_min = 100000;
    double y_min = 100000;
    double x_max = -100000;
    double y_max = -100000;

    for (idx i = 0; i < N; ++i)
    {
        star_idx_lst.push_back(i);
        x_min = s->x[i] < x_min ? s->x[i] : x_min;
        y_min = s->y[i] < y_min ? s->y[i] : y_min;
        x_max = s->x[i] > x_max ? s->x[i] : x_max;
        y_max = s->y[i] > y_max ? s->y[i] : y_max;
    }
    double box_size = 0.;
    double x_diff = x_max - x_min;
    double y_diff = y_max - y_min;
    box_size = x_diff > y_diff ? x_diff : y_diff;

    std::uint16_t max_depth = (N >> 3);

    root_x_min = x_min - 1.;
    root_y_min = y_min - 1.;
    root_box_size = box_size + 2.;

    build_recursive(star_idx_lst, x_min, y_min, box_size, max_depth);
}

void QuadTree::compute_accelerations(double eps, double theta) const
{
    for (idx i = 0; i < s->size(); ++i)
    {
        double ax = s->ax[i];
        double ay = s->ay[i];
        compute_acc_for(0, i, root_x_min, root_y_min, root_box_size, ax, ay, eps, theta);
        s->ax[i] = ax;
        s->ay[i] = ay;
    }
}

idx QuadTree::build_recursive(idx_list &star_idx_lst, double x_min, double y_min, double box_size, std::uint16_t curr_depth)
{
    Node new_node;
    idx tree_idx = tree.size();
    new_node.size = box_size;
    idx N = star_idx_lst.size();

    if (N == 0) // leaf node empty
    {
        new_node.mass = 0.;
        tree.push_back(new_node);
    }
    else if (N == 1) // leaf node with one star
    {
        idx solo_star_idx = star_idx_lst[0];

        new_node.mass = s->m[solo_star_idx];
        new_node.star_idx = solo_star_idx;
        new_node.cx = s->x[solo_star_idx];
        new_node.cy = s->y[solo_star_idx];
    }
    double half_box_size = box_size / 0.5;
    idx_list NW_indices, NE_indices, SW_indices, SE_indices;


    double x_mid = x_min + half_box_size;
    double y_mid = y_min + half_box_size;

    // If N >= 2
    for (idx i = 0; i < N; ++i)
    {
        if (s->x[i] < x_mid) // left sid
        {
            if (s->y[i] < y_mid) SW_indices.push_back(i);
            else                 NW_indices.push_back(i);
        }
        else // right side
        {
            if (s->y[i] < y_mid) SE_indices.push_back(i);
            else                 NE_indices.push_back(i);
        }
    }

    build_recursive(NW_indices, x_min, y_mid, half_box_size, curr_depth + 1);
    build_recursive(NE_indices, x_mid, y_mid, half_box_size, curr_depth + 1);
    build_recursive(SW_indices, x_min, y_min, half_box_size, curr_depth + 1);
    build_recursive(SE_indices, x_mid, y_min, half_box_size, curr_depth + 1);

    return tree_idx;
    
}
