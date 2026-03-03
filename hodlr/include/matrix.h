#ifndef MATRIX_H
#define MATRIX_H

#include <vector>


class Random;

class Matrix 
{
    public:
        std::vector<double> mat;

        Matrix(std::size_t, std::size_t);
        std::size_t getSize() const;
        double operator() (const std::size_t, const std::size_t) const;

        std::size_t getNumRows() const;
        std::size_t getNumCols() const;

        void fillRandom(Random&);

    private:
        std::size_t num_rows;
        std::size_t num_cols;
};

#endif //!MATRIX_H
