#include <iostream>
#include <cmath>
#include "../include/hodlr.h"

void printMetaData(HODLR_Matrix&);

int main(int argc, char** argv)
{
   int n = std::stoi(argv[1]); // Size of leaf nodes: n x n
   int k = std::stoi(argv[2]); // Rank of low-rank off-diagonal blocks
   int L = std::stoi(argv[3]); // Number of levels
   unsigned int seed = std::stoi(argv[4]); // random seed

   Random rng(seed);
   HODLR_Matrix hodlr(n, k, L);
   hodlr.fillWithRandomData(rng);
   hodlr.initializeTree();

   hodlr.printMetaData();
   hodlr.printTreeData();


    return 0;
}

