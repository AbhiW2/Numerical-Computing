
#ifndef GAUSSIAN_ELIMINATION_HPP
#define GAUSSIAN_ELIMINATION_HPP

#include "LinearSystemSolver.hpp"
#include <vector>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

class GaussianElimination : public LinearSystemSolver {
private:
    static constexpr double EPS = 1e-10;

    void forwardElimination();
    std::vector<double> backSubstitution();

public:
    explicit GaussianElimination(int n);

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    std::vector<double> solve() override;
};

#endif