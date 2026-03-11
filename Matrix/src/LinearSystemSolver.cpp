#include "../Include/LinearSystemSolver.hpp"

LinearSystemSolver::LinearSystemSolver(int n)
    : Matrix(n, n + 1)
{}

LinearSystemSolver::LinearSystemSolver(int r, int c)
    : Matrix(r, c)
{}