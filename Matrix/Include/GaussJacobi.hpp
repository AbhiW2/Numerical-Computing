#ifndef GAUSS_JACOBI_HPP
#define GAUSS_JACOBI_HPP

#include "LinearSystemSolver.hpp"
#include <vector>
#include <cmath>

class GaussJacobi : public LinearSystemSolver {
private:
    int    maxIter;
    double tol;

public:
    GaussJacobi(int n, int maxIter = 1000, double tol = 1e-6);

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    std::vector<double> solve() override;
};

#endif