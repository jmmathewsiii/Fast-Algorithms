#include <iostream>
#include "../include/hodlr.h"

HODLR_Matrix::HODLR_Matrix(std::size_t leaf_size, std::size_t low_rank, std::size_t num_levels)
    : n(leaf_size), k(low_rank), L(num_levels), N((static_cast<std::size_t>(1) << num_levels) * leaf_size)
{
    leaf_vals.reserve(n * N);
    u_vals.reserve(L * N * k);
    vt_vals.reserve(L * N * k);
    tree.reserve((static_cast<std::size_t>(1) << (num_levels + 1)) - 1);
}

std::size_t HODLR_Matrix::getNumRows() const { return N; }

std::size_t HODLR_Matrix::getLowRank() const { return k; }

std::size_t HODLR_Matrix::getLeafSize() const { return n; }

std::size_t HODLR_Matrix::getNumLevels() const { return L; }

std::size_t HODLR_Matrix::getUPoolSize() const { return u_vals.size(); }
std::size_t HODLR_Matrix::getVTPoolSize() const { return vt_vals.size(); }
std::size_t HODLR_Matrix::getLeafPoolSize() const { return leaf_vals.size(); }
std::size_t HODLR_Matrix::getTreeSize() const { return tree.size(); }

void HODLR_Matrix::fillWithRandomData(Random &rng) {
    for (std::size_t i = 0; i < n * N; ++i)
        leaf_vals.push_back(rng.uniform());
    for (std::size_t i = 0; i < L * N * k; ++i) {
        u_vals.push_back(rng.uniform());
        vt_vals.push_back(rng.uniform());
    }
}

void HODLR_Matrix::initializeTree() {
    for (std::size_t l = 0; l < L; ++l)
    {
        const std::size_t level_start_idx = (static_cast<std::size_t>(1) << l) - 1;
        const std::size_t level_end_idx = (static_cast<std::size_t>(1) << (l + 1)) - 1;
        const std::size_t mat_block_size = N / (static_cast<std::size_t>(1) << l);

        const std::size_t low_rank_mat_size = (N * k) / (static_cast<std::size_t>(1) << (l + 1));
        const std::size_t Nkl = N * k * l;

        for (std::size_t idx = level_start_idx; idx < level_end_idx; ++idx) {
            const std::size_t it = idx - level_start_idx; // Same as idx but ignoring the offset
            Node node;
            node.level = l;
            node.start = it * mat_block_size;
            node.end = node.start + mat_block_size;
            node.isLeaf = false;

            node.left_child = (idx * 2) + 1;
            node.right_child = node.left_child + 1;
            node.U1 = Nkl + (2 * it * low_rank_mat_size);
            node.U2 = node.U1 + low_rank_mat_size;

            node.VT1 = Nkl + (2 * it * low_rank_mat_size);
            node.VT2 = node.VT1 + low_rank_mat_size;

            tree.push_back(node);
        }
    }
    const std::size_t leaf_end = (static_cast<std::size_t>(1) << L);

    for (std::size_t i = 0; i < leaf_end; ++i)
    {
        const std::size_t n_sq = n * n;

        Node node;
        node.level = L;
        node.start = i * n;
        node.end = node.start + n;
        node.isLeaf = true;
        node.A = i * n_sq;

        tree.push_back(node);
    }
}

void HODLR_Matrix::printMetaData() {
   std::cout << "\n        Low rank: " << getLowRank()
             << "\nLeaf matrix size: " << getLeafSize()
             << "\n      Leaf level: " << getNumLevels()
             << "\nFull matrix size: " << getNumRows() 
             << "\n     U Pool size: " << getUPoolSize()
             << "\n    VT Pool size: " << getVTPoolSize()
             << "\n  Leaf Pool size: " << getLeafPoolSize()
             << "\n       Tree size: " << getTreeSize() 
             << "\n\n";
}

void HODLR_Matrix::printTreeData() {
    for (std::size_t i = 0; i < tree.size(); ++i) {
        Node c = tree.at(i);

        std::cout << "\n             Tree Idx: " << i
                  << "\n                Level: " << c.level
                  << "\n               isLeaf: " << c.isLeaf
                  << "\nFull Matrix Range: [" << c.start << ", " << c.end << ")";
        if (!c.isLeaf) {
            std::cout << "\n       Left Child Idx: " << c.left_child
                      << "\n      Right Child Idx: " << c.right_child
                      << "\n   U12 Pool Start Idx: " << c.U1
                      << "\n   U21 Pool Start Idx: " << c.U2
                      << "\n  VT12 Pool Start Idx: " << c.VT1
                      << "\n  VT21 Pool Start Idx: " << c.VT2
                      << "\n";
        }
        else {
            std::cout <<"\n  Leaf Pool Start Idx: " << c.A << "\n";
        }
    }
}
