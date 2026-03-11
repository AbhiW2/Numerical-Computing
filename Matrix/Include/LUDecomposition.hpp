#ifndef LUDECOMPOSITION_HPP
#define LUDECOMPOSITION_HPP

#include "LinearSystemSolver.hpp"
#include <vector>

class LUDecomposition : public LinearSystemSolver {
private:
    Matrix L;
    Matrix U;

public:
    LUDecomposition(int n);

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    void decompose();
    std::vector<double> forwardSub(const std::vector<double>& b);
    std::vector<double> backSub(const std::vector<double>& y);
    std::vector<double> solve() override;
};

#endif // LUDECOMPOSITION_HPP