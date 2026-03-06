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

VD HODLR_Matrix::createFullMatrix() {
    VD full_mat(N * N);
    std::size_t leaf_start_idx = static_cast<std::size_t>(1) << L;
    for (std::size_t tree_idx = 0; tree_idx < leaf_start_idx; ++tree_idx) {
        Node c = tree[tree_idx];
        VD UVT1(c.uv_length * c.uv_length);
        VD UVT2(c.uv_length * c.uv_length);
        for (std::size_t i = 0; i < c.uv_length; ++i) {
            for (std::size_t j = 0; j < c.uv_length; ++j) {
                for (std::size_t r = 0; r < k; ++r) {
                    UVT1[(i * c.uv_length) + j] += 
                        u_vals[c.U1 + (i * k) + r] * vt_vals[c.VT1 + (r * c.uv_length) + j];
                    UVT2[(i * c.uv_length) + j] += 
                        u_vals[c.U2 + (i * k) + r] * vt_vals[c.VT2 + (r * c.uv_length) + j];
                }
            }
        }
        for (std::size_t i = 0; i < c.uv_length; ++i) {
            for (std::size_t j = 0; j < c.uv_length; ++j) {
                full_mat[(c.start * N) + c.start + c.uv_length + (i * N) + j] = UVT1[(i * c.uv_length) + j];
                full_mat[(c.start * N) + c.start + (c.uv_length * N) + (i * N) + j] = UVT2[(i * c.uv_length) + j];
            }
        }
    }
    for (std::size_t tree_idx = leaf_start_idx; tree_idx < tree.size(); ++tree_idx) {
        Node c = tree[tree_idx];
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                full_mat[(c.start * N) + c.start + ((i * N) + j)] = leaf_vals[c.A + ((i * n) + j)];
    }

    return full_mat;
}

VD HODLR_Matrix::MatVec(VD &x) {
    if (x.size() != N) {
        std::cerr << "HODLR Matrix Size: " << N << " x " << N
                  << "\nVector size: " << x.size() << "\n";
        VD fail(0);
        return fail;
    }

    VD b(N);
    std::size_t tree_size = tree.size();
    std::size_t leaf_start_idx = static_cast<std::size_t>(1) << L;

    for (std::size_t tree_idx = 0; tree_idx < leaf_start_idx; ++tree_idx) {
       Node c = tree.at(tree_idx); 

       VD tmp1(k);                                       // Intermediate vector storing VT12 * x
       VD tmp2(k);                                       // Intermediate vector storing VT21 * x
       std::size_t x1_start_idx = c.start + c.uv_length; // Start index into x for multiplication by VT12
       std::size_t x2_start_idx = c.start;               // Start index into x for multiplication by VT21
       std::size_t b1_start_idx = c.start;               // Start index into b for multiplication by U12tmp
       std::size_t b2_start_idx = c.start + c.uv_length; // Start index into b for multiplication by U21tmp

       for (std::size_t i = 0; i < k; ++i) {
           for (std::size_t j = 0; j < c.uv_length; ++j) {
               tmp1[i] += vt_vals[c.VT1 + ((i * c.uv_length) + j)] * x[x1_start_idx + j];
               tmp2[i] += vt_vals[c.VT2 + ((i * c.uv_length) + j)] * x[x2_start_idx + j];
           }
       }
       for (std::size_t i = 0; i < c.uv_length; ++i) {
           for (std::size_t j = 0; j < k; ++j) {
               b[b1_start_idx + i] += u_vals[c.U1 + ((i * k) + j)] * tmp1[j];
               b[b2_start_idx + i] += u_vals[c.U2 + ((i * k) + j)] * tmp2[j];
           }
       }
    }

    for (std::size_t tree_idx = leaf_start_idx; tree_idx < tree_size; ++tree_idx) {
        Node c = tree.at(tree_idx);
        for (std::size_t i = 0; i < n; ++i) 
            for (std::size_t j = 0; j < n; ++j)
                b[c.start + i] += leaf_vals[c.A + ((i * n) + j)] * x[c.start + j];
    }

    return b;
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
            node.length = mat_block_size;
            node.uv_length = mat_block_size / static_cast<std::size_t>(2);
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
