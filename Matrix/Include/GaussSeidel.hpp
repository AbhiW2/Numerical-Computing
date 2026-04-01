#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

#include "LinearSystemSolver.hpp"
#include <vector>
#include <cmath>

class GaussSeidel : public LinearSystemSolver {
private:
    int    maxIter;
    double tol;

public:
    GaussSeidel(int n, int maxIter = 1000, double tol = 1e-6);

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    std::vector<double> solve() override;
};

#endif