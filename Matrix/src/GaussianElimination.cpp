#include "../Include/GaussianElimination.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace std;

GaussianElimination::GaussianElimination(int n)
    : LinearSystemSolver(n, n + 1)
{}

void GaussianElimination::loadAugmentedMatrix(ifstream& Afile, ifstream& /*bfile*/) {
    loadFromAugmentedFile(Afile);
}

void GaussianElimination::forwardElimination() {
    int n = rows;
    for (int k = 0; k < n; k++) {
        // Partial pivoting
        int maxRow = k;
        double maxVal = fabs((*this)(k, k));
        for (int i = k + 1; i < n; i++) {
            if (fabs((*this)(i, k)) > maxVal) {
                maxVal = fabs((*this)(i, k));
                maxRow = i;
            }
        }
        if (maxRow != k)
            for (int j = 0; j <= n; j++)
                swap((*this)(k, j), (*this)(maxRow, j));

        if (fabs((*this)(k, k)) < EPS) continue;

        for (int i = k + 1; i < n; i++) {
            double factor = (*this)(i, k) / (*this)(k, k);
            for (int j = k; j <= n; j++)
                (*this)(i, j) -= factor * (*this)(k, j);
        }
    }
}

vector<double> GaussianElimination::backSubstitution() {
    int n = rows;
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        if (fabs((*this)(i, i)) < EPS) {
            cout << "Warning: Free variable at row " << i << " (set to 0)\n";
            continue;
        }
        double sum = (*this)(i, n);
        for (int j = i + 1; j < n; j++)
            sum -= (*this)(i, j) * x[j];
        x[i] = sum / (*this)(i, i);
    }
    return x;
}

vector<double> GaussianElimination::solve() {
    if (rows <= 0) throw runtime_error("Invalid matrix size");
    forwardElimination();
    return backSubstitution();
}