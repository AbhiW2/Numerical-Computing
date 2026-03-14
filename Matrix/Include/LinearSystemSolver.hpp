// #ifndef LINEAR_SYSTEM_SOLVER_HPP
// #define LINEAR_SYSTEM_SOLVER_HPP

// #include "Matrix.hpp"
// #include <vector>
// #include <fstream>

// class LinearSystemSolver : public Matrix {
// public:
//     LinearSystemSolver(int n);        // n equations → augmented matrix n×(n+1)
//     LinearSystemSolver(int r, int c); // general constructor
//     virtual ~LinearSystemSolver() = default;

//     virtual void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) = 0;
//     virtual std::vector<double> solve() = 0;
// };

// #endif // LINEAR_SYSTEM_SOLVER_HPP








#ifndef LINEAR_SYSTEM_SOLVER_HPP
#define LINEAR_SYSTEM_SOLVER_HPP

#include "Matrix.hpp"
#include <vector>
#include <fstream>

class LinearSystemSolver : public Matrix {
public:
    LinearSystemSolver(int n);        // n equations → augmented matrix n×(n+1)
    LinearSystemSolver(int r, int c); // general constructor
    virtual ~LinearSystemSolver() = default;

    virtual void loadAugmentedMatrix(std::ifstream& Afile, std::ifstream& bfile) = 0;
    virtual std::vector<double> solve() = 0;
};

#endif // LINEAR_SYSTEM_SOLVER_HPP