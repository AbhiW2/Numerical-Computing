#ifndef GAUSSIAN_ELIMINATION_HPP
#define GAUSSIAN_ELIMINATION_HPP

#include "Matrix.hpp"
#include <vector>
#include <fstream>

class GaussianElimination : public Matrix {
public:
    // Constructor
    GaussianElimination(int r, int c);

    // Combine A (49l.txt) and b (49r.txt) → augmented matrix
    void loadAugmentedMatrix(std::ifstream &Afile, std::ifstream &bfile);

    // Gaussian elimination steps
    void forwardElimination();
    std::vector<double> backSubstitution();

    // Solve system
    std::vector<double> solve();
};

#endif