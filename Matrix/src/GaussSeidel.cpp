#include "../Include/GaussSeidel.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

GaussSeidel::GaussSeidel(int n, int maxIter, double tol)
    : LinearSystemSolver(n), maxIter(maxIter), tol(tol) {}

void GaussSeidel::loadAugmentedMatrix(ifstream& Afile, ifstream& /*bfile*/) {
    loadFromAugmentedFile(Afile);
}


// ── solve ─────────────────────────────────────────────────────────────────────
vector<double> GaussSeidel::solve() {
    int n = rows;
    int augC = cols;

    // Step 1: Diagonal dominance
    if (!checkDiagDominant()) {
        cout << "Gauss-Seidel: matrix is NOT diagonally dominant.\n";
        cout << "Attempting reordering...\n";

        makeDiagDominant();

        if (checkDiagDominant())
            cout << "Reordering successful\n\n";
        else
            cout << "Still not diagonally dominant\n\n";
    }

    // Step 2: identify free variables
    vector<bool> free_var(n, false);
    int free_count = 0;

    for (int i = 0; i < n; i++) {
        if (fabs((*this)(i, i)) < 1e-12) {
            free_var[i] = true;
            free_count++;
        }
    }

    if (free_count > 0)
        cout << "Note: " << free_count << " free variables\n\n";

    // Step 3: iteration
    vector<double> x(n, 0.0);

    for (int iter = 0; iter < maxIter; iter++) {
        double error = 0.0;

        for (int i = 0; i < n; i++) {
            if (free_var[i]) continue;

            double diag = (*this)(i, i);
            double sum  = (*this)(i, n);

            for (int j = 0; j < n; j++) {
                if (j != i)
                    sum -= (*this)(i, j) * x[j];
            }

            double newVal = sum / diag;

            if (isnan(newVal) || isinf(newVal))
                throw runtime_error("Divergence detected");

            error += fabs(newVal - x[i]);
            x[i] = newVal;
        }

        if (error < tol) {
            cout << "Converged in " << iter + 1 << " iterations\n";
            return x;
        }
    }

    cout << "Max iterations reached\n";
    return x;
}