#ifndef GAUSSIAN_ELIMINATION_HPP
#define GAUSSIAN_ELIMINATION_HPP

#include "LinearSystemSolver.hpp"
#include <vector>
#include <fstream>
#include <cmath>

class GaussianElimination : public LinearSystemSolver {
private:
    static constexpr double EPS = 1e-12;

public:
    GaussianElimination(int n);

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    void forwardElimination();
    std::vector<double> backSubstitution();
    std::vector<double> solve() override;
};

#endif