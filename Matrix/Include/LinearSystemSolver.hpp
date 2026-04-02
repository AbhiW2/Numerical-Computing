
#ifndef LINEAR_SYSTEM_SOLVER_HPP
#define LINEAR_SYSTEM_SOLVER_HPP

#include "Matrix.hpp"
#include <vector>
#include <fstream>

class LinearSystemSolver : public Matrix {
public:
    LinearSystemSolver(int n);
    LinearSystemSolver(int r, int c);
    virtual ~LinearSystemSolver() = default;

    virtual void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) = 0;
    virtual std::vector<double> solve() = 0;

protected:
    // Reads "n n+1" header then fills the n×(n+1) augmented matrix [A|b]
    void loadFromAugmentedFile(std::ifstream& Afile);
    bool checkDiagDominant();
    void makeDiagDominant();
};

#endif