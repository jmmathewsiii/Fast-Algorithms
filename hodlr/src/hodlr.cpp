#include "../include/random.h"
#include <cstddef>
#include "../include/hodlr.h"

HODLR_Matrix::HODLR_Matrix(std::size_t leaf_size, std::size_t low_rank, std::size_t num_levels)
    : n(leaf_size), k(low_rank), L(num_levels), N((static_cast<std::size_t>(1) << num_levels) * leaf_size)
{
    leaf_vals.reserve(n * N);
    u_vals.reserve(L * N * k);
    vt_vals.reserve(L * N * k);
    tree.reserve((static_cast<std::size_t>(1) << (num_levels + 1)) - 1);
}

std::size_t HODLR_Matrix::getNumRows() const {
    return N;
}

std::size_t HODLR_Matrix::getLeafSize() const {
    return n;
}

std::size_t HODLR_Matrix::getNumLevels() const {
    return L;
}

void HODLR_Matrix::fillWithRandomData(Random &rng) {
    for (std::size_t i = 0; i < leaf_vals.capacity(); ++i)
        leaf_vals.push_back(rng.uniform());
    for (std::size_t i = 0; i < u_vals.capacity(); ++i) {
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
        const std::size_t half_block_size = mat_block_size / 2;

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
            node.U1 = Nkl + (it * low_rank_mat_size);
            node.U2 = node.U1 + low_rank_mat_size;

            node.VT1 = Nkl + (it * low_rank_mat_size);
            node.VT2 = node.U1 + low_rank_mat_size;

            tree.push_back(node);
        }
    }
    const std::size_t leaf_end = (1 << L);

    for (std::size_t i = 0; i < leaf_end; ++i)
    {
        const std::size_t n_sq = n * n;

        Node node;
        node.level = L;
        node.start = i * n_sq;
        node.end = node.start + n_sq;
        node.isLeaf = true;
        node.A = i * n_sq;

        tree.push_back(node);
    }
}
