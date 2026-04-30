#include "../include/quadtree_cpu.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>

QuadTree::QuadTree(Stars *s_) : s(s_) {}

const Tree& QuadTree::getTree() const { return tree; }

void QuadTree::build()
{
    tree.clear();
    tree.reserve(s->size() << 1);
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

    double boundary_buffer = box_size * 0.1;

    root_x_min = x_min - boundary_buffer;
    root_y_min = y_min - boundary_buffer;
    root_box_size = box_size + 2 * boundary_buffer;

    Node root;
    root.size = root_box_size;
    tree.push_back(root);

    build_recursive(0, star_idx_lst, root_x_min, root_y_min, root_box_size, 0);
}


void QuadTree::build_recursive(idx curr_idx, idx_list &star_idx_lst, double x_min, double y_min, double box_size, std::uint16_t curr_depth)
{
    idx N = star_idx_lst.size();

    if (N == 0) // leaf node empty
    {
        tree[curr_idx].mass = 0.;
        tree[curr_idx].cx = 0.;
        tree[curr_idx].cy = 0.;
        tree[curr_idx].first_child = NO_CHILD;
        tree[curr_idx].star_idx = NO_CHILD;
    }
    else if (N == 1) // leaf node with one star
    {
        idx solo_star_idx = star_idx_lst[0];

        tree[curr_idx].mass = s->m[solo_star_idx];
        tree[curr_idx].star_idx = solo_star_idx;
        tree[curr_idx].cx = s->x[solo_star_idx];
        tree[curr_idx].cy = s->y[solo_star_idx];
        tree[curr_idx].first_child = NO_CHILD;
    }
    else // non-leaf node
    {
        if (curr_depth >= MAX_DEPTH) throw std::runtime_error("Exceeded max depth");
        double half_box_size = box_size * 0.5;
        idx_list NW_indices, NE_indices, SW_indices, SE_indices;

        double x_mid = x_min + half_box_size;
        double y_mid = y_min + half_box_size;

        // If N >= 2
        for (idx i = 0; i < N; ++i)
        {
            idx p = star_idx_lst[i];
            if (s->x[p] < x_mid) // left side
            {
                if (s->y[p] < y_mid) SW_indices.push_back(p);
                else                 NW_indices.push_back(p);
            }
            else // right side
            {
                if (s->y[p] < y_mid) SE_indices.push_back(p);
                else                 NE_indices.push_back(p);
            }
        }

        idx first_child = tree.size();
        tree.resize(first_child + 4);

        build_recursive(first_child    , NW_indices, x_min, y_mid, half_box_size, curr_depth + 1);
        build_recursive(first_child + 1, NE_indices, x_mid, y_mid, half_box_size, curr_depth + 1);
        build_recursive(first_child + 2, SW_indices, x_min, y_min, half_box_size, curr_depth + 1);
        build_recursive(first_child + 3, SE_indices, x_mid, y_min, half_box_size, curr_depth + 1);

        double total_mass = tree[first_child    ].mass + tree[first_child + 1].mass 
                          + tree[first_child + 2].mass + tree[first_child + 3].mass;

        double weighted_cx = tree[first_child    ].cx * tree[first_child    ].mass
                           + tree[first_child + 1].cx * tree[first_child + 1].mass
                           + tree[first_child + 2].cx * tree[first_child + 2].mass
                           + tree[first_child + 3].cx * tree[first_child + 3].mass;

        double weighted_cy = tree[first_child    ].cy * tree[first_child    ].mass
                           + tree[first_child + 1].cy * tree[first_child + 1].mass
                           + tree[first_child + 2].cy * tree[first_child + 2].mass
                           + tree[first_child + 3].cy * tree[first_child + 3].mass;

       tree[curr_idx].mass = total_mass;
       tree[curr_idx].cx = weighted_cx / total_mass;
       tree[curr_idx].cy = weighted_cy / total_mass;
       tree[curr_idx].star_idx = NO_CHILD;
       tree[curr_idx].first_child = first_child;
    }
}

void QuadTree::compute_accelerations(double eps, double theta) const
{
    for (idx i = 0; i < s->size(); ++i)
    {
        double ax = 0.;
        double ay = 0.;
        compute_acc_for(0, i, ax, ay, eps, theta);
        s->ax[i] = ax;
        s->ay[i] = ay;
    }
}

void QuadTree::compute_acc_for(idx node_idx, idx star_idx, double &ax, double &ay, double eps, double theta) const
{
    if (tree[node_idx].mass == 0.) return; // empty leaf
    if (tree[node_idx].star_idx == star_idx) return; // self-interaction

    double my_x = s->x[star_idx];
    double my_y = s->y[star_idx];
    double dx = my_x - tree[node_idx].cx;
    double dy = my_y - tree[node_idx].cy;

    double r_sq = dx * dx + dy * dy + (eps * eps);
    double distance = std::sqrt(r_sq);
        
    if (tree[node_idx].first_child == NO_CHILD) // single particle
    {
        double r3_inv = 1. / (r_sq * distance);
        ax -= G * dx * tree[node_idx].mass * r3_inv;
        ay -= G * dy * tree[node_idx].mass * r3_inv;
    }
    else if (tree[node_idx].size / distance < theta) // node within tolerance
    {
        double r3_inv = 1. / (r_sq * distance);
        ax -= G * dx * tree[node_idx].mass * r3_inv;
        ay -= G * dy * tree[node_idx].mass * r3_inv;
    }
    else // node outside tolerance
    {
        idx first_child = tree[node_idx].first_child;
        compute_acc_for(first_child, star_idx, ax, ay, eps, theta);
        compute_acc_for(first_child + 1, star_idx, ax, ay, eps, theta);
        compute_acc_for(first_child + 2, star_idx, ax, ay, eps, theta);
        compute_acc_for(first_child + 3, star_idx, ax, ay, eps, theta);
    }
}

void QuadTree::write_divisions(std::ostream& os) const
{
    if (tree.empty()) return;

    double x_max = root_x_min + root_box_size;
    double y_max = root_y_min + root_box_size;
    os << root_x_min << " " << root_y_min << "\n";
    os << x_max      << " " << root_y_min << "\n";
    os << x_max      << " " << y_max      << "\n";
    os << root_x_min << " " << y_max      << "\n";
    os << root_x_min << " " << root_y_min << "\n";
    os << "\n";

    write_divisions_recursive(0, root_x_min, root_y_min, root_box_size, os);
}

void QuadTree::write_divisions_recursive(idx node_idx, double x_min, double y_min, double box_size, std::ostream& os) const
{
    if (tree[node_idx].first_child == NO_CHILD) return;

    double half = box_size * 0.5;
    double x_mid = x_min + half;
    double y_mid = y_min + half;
    double x_max = x_min + box_size;
    double y_max = y_min + box_size;

    os << x_mid << " " << y_min << "\n";
    os << x_mid << " " << y_max << "\n";
    os << "\n";

    os << x_min << " " << y_mid << "\n";
    os << x_max << " " << y_mid << "\n";
    os << "\n";

    idx fc = tree[node_idx].first_child;
    write_divisions_recursive(fc    , x_min, y_mid, half, os);  // NW
    write_divisions_recursive(fc + 1, x_mid, y_mid, half, os);  // NE
    write_divisions_recursive(fc + 2, x_min, y_min, half, os);  // SW
    write_divisions_recursive(fc + 3, x_mid, y_min, half, os);  // SE
}














