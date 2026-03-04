#ifndef HODLR_H
#define HODLR_H

#include <vector>
#include "node.h"

using VD = std::vector<double>;

class Random;

class HODLR_Matrix {
    private:
        const std::size_t n; // number of rows/columns of a leaf node
        const std::size_t k; // rank of the off-diagnonal entries
        const std::size_t L; // Number of levels in the HODLR matrix
        const std::size_t N; // Full size of the matrix

        VD leaf_vals;
        VD u_vals;
        VD vt_vals; // V transpose (cache-friendly)

        std::vector<Node> tree;

    public:
        HODLR_Matrix(std::size_t, std::size_t, std::size_t);
        std::size_t getNumRows() const;
        std::size_t getLowRank() const;
        std::size_t getLeafSize() const;
        std::size_t getNumLevels() const;

        void fillWithRandomData(Random&);

        void initializeTree();
};

#endif // !HODLR_H
