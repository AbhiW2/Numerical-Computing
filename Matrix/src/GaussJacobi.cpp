#include "../Include/GaussJacobi.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

GaussJacobi::GaussJacobi(int n, int maxIter, double tol)
    : LinearSystemSolver(n), maxIter(maxIter), tol(tol) {}

void GaussJacobi::loadAugmentedMatrix(ifstream& Afile, ifstream& /*bfile*/) {
    loadFromAugmentedFile(Afile);
}

// ── check diagonal dominance (A part only, ignore b column) ──────────────────
static bool checkDiagDominant(double** data, int n) {
    for (int i = 0; i < n; i++) {
        double diag = fabs(data[i][i]);
        double sum  = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) sum += fabs(data[i][j]);
        if (diag <= sum) return false;
    }
    return true;
}

// ── greedy row-swap: for each column k, assign unassigned row with max |A[i][k]| ──
static void makeDiagDominant(double** data, int n, int augCols) {
    vector<bool> used(n, false);
    vector<int>  order(n, -1);

    for (int k = 0; k < n; k++) {
        int    best    = -1;
        double bestVal = -1.0;
        for (int i = 0; i < n; i++) {
            if (used[i]) continue;
            double v = fabs(data[i][k]);
            if (v > bestVal) { bestVal = v; best = i; }
        }
        if (best == -1)
            throw runtime_error("Diagonal dominance reordering: no valid pivot found");
        order[k] = best;
        used[best] = true;
    }

    double** tmp = new double*[n];
    for (int k = 0; k < n; k++) {
        tmp[k] = new double[augCols];
        for (int j = 0; j < augCols; j++)
            tmp[k][j] = data[order[k]][j];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < augCols; j++)
            data[i][j] = tmp[i][j];
        delete[] tmp[i];
    }
    delete[] tmp;
}

// ── solve ─────────────────────────────────────────────────────────────────────
vector<double> GaussJacobi::solve() {
    int n    = rows;
    int augC = cols;   // n+1

    // Step 1: diagonal dominance check and reorder
    if (!checkDiagDominant(data, n)) {
        cout << "Gauss-Jacobi: matrix is NOT diagonally dominant.\n"
             << "              Attempting greedy row reordering...\n";
        makeDiagDominant(data, n, augC);

        if (checkDiagDominant(data, n))
            cout << "              Reordering successful - now diagonally dominant.\n\n";
        else
            cout << "              Warning: diagonal dominance could NOT be achieved.\n"
                 << "              Continuing anyway (results may vary).\n\n";
    } else {
        cout << "Gauss-Jacobi: matrix is diagonally dominant. No reordering needed.\n\n";
    }

    // Step 2: identify free variables (zero diagonal rows → x[i] = 0, skipped)
    vector<bool> free_var(n, false);
    int free_count = 0;
    for (int i = 0; i < n; i++) {
        if (fabs((*this)(i, i)) < 1e-12) {
            free_var[i] = true;
            free_count++;
        }
    }
    if (free_count > 0)
        cout << "Note: " << free_count << " zero-diagonal row(s) detected — "
             << "treating as free variables (set to 0).\n\n";

    // Step 3: Jacobi iteration — uses OLD x values for all updates
    vector<double> x(n, 0.0), x_new(n, 0.0);

    for (int iter = 0; iter < maxIter; iter++) {
        double error = 0.0;

        for (int i = 0; i < n; i++) {
            if (free_var[i]) {         // skip free variables
                x_new[i] = 0.0;
                continue;
            }

            double diag = (*this)(i, i);
            double sum  = (*this)(i, n);   // b[i]
            for (int j = 0; j < n; j++)
                if (j != i) sum -= (*this)(i, j) * x[j];   // OLD x

            x_new[i] = sum / diag;

            if (isnan(x_new[i]) || isinf(x_new[i]))
                throw runtime_error(
                    "Gauss-Jacobi diverged (NaN/Inf) at iteration " + to_string(iter + 1));

            error += fabs(x_new[i] - x[i]);
        }

        x = x_new;   // update all at once

        if (error < tol) {
            cout << "Gauss-Jacobi converged in " << iter + 1 << " iterations.\n";
            return x;
        }
    }

    cout << "Gauss-Jacobi reached max iterations (" << maxIter
         << ") — returning best approximation.\n";
    return x;
}