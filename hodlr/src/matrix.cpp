#include "../include/matrix.h"
#include "../include/random.h"



Matrix::Matrix(std::size_t n_rows, std::size_t n_cols)
    : num_rows(n_rows), num_cols(n_cols)   
{
    mat.reserve(n_rows * n_cols);
}

std::size_t Matrix::getSize() const {
    return mat.size();
}

double Matrix::operator() (const std::size_t i, const std::size_t j) const {
    return mat.at(i * num_rows + j);
}

std::size_t Matrix::getNumRows() const {
    return num_rows;
}

void Matrix::fillRandom(Random &rng) {
    for (std::size_t i = 0; i < mat.size(); ++i){
        mat[i] = rng.uniform();
    }
}
