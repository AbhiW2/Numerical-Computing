#ifndef LU_DECOMPOSITION_HPP
#define LU_DECOMPOSITION_HPP

#include "LinearSystemSolver.hpp"
#include <vector>
#include <fstream>
#include <cmath>
#include <stdexcept>

enum class LUMethod {
    DOOLITTLE,
    CROUT,
    CHOLESKY
};

class LUDecomposition : public LinearSystemSolver {
public:
    explicit LUDecomposition(int n, LUMethod m = LUMethod::DOOLITTLE);
    ~LUDecomposition() override = default;

    void     setMethod(LUMethod newMethod);
    LUMethod getMethod() const { return method; }

    void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) override;
    std::vector<double> solve() override;

private:
    static constexpr double EPS = 1e-10;

    Matrix   L;
    Matrix   U;
    LUMethod method;

    std::vector<int> doolittle();
    std::vector<int> crout();
    std::vector<int> cholesky();
    std::vector<int> decompose();

    std::vector<double> forwardSubstitution(const Matrix& Lmat,
                                            const std::vector<double>& b);
    std::vector<double> backSubstitution   (const Matrix& Umat,
                                            const std::vector<double>& y);
};

#endif