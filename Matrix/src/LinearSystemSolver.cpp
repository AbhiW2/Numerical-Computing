#include "../Include/LinearSystemSolver.hpp"
#include <stdexcept>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

LinearSystemSolver::LinearSystemSolver(int n)
    : Matrix(n, n + 1)
{}

LinearSystemSolver::LinearSystemSolver(int r, int c)
    : Matrix(r, c)
{}

// Reads the 49l.txt augmented-matrix file.
// Format: first line = "n  n+1", then n rows of (n+1) values.
// Last value in each row is b[i]. No separate bfile needed.
void LinearSystemSolver::loadFromAugmentedFile(ifstream& Afile)
{
    int file_rows, file_cols;
    if (!(Afile >> file_rows >> file_cols))
        throw runtime_error("Cannot read matrix dimensions from file");

    if (file_rows != rows || file_cols != cols)
        throw runtime_error(
            "Size mismatch: solver expects " +
            to_string(rows) + "x" + to_string(cols) +
            " but file has " + to_string(file_rows) +
            "x" + to_string(file_cols));

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (!(Afile >> data[i][j]))
                throw runtime_error(
                    "Error reading element (" + to_string(i) +
                    "," + to_string(j) + ") from file");

    cout << "Augmented matrix [A|b] loaded ("
         << rows << " equations, " << rows << " unknowns).\n\n";
}


bool LinearSystemSolver::checkDiagDominant() {
    for (int i = 0; i < rows; i++) {
        double diag = fabs(data[i][i]);
        double sum  = 0.0;

        for (int j = 0; j < rows; j++) {
            if (j != i) sum += fabs(data[i][j]);
        }

        if (diag <= sum) return false;
    }
    return true;
}

void LinearSystemSolver::makeDiagDominant() {
    vector<bool> used(rows, false);
    vector<int> order(rows, -1);

    for (int k = 0; k < rows; k++) {
        int best = -1;
        double bestVal = -1.0;

        for (int i = 0; i < rows; i++) {
            if (used[i]) continue;

            double v = fabs(data[i][k]);
            if (v > bestVal) {
                bestVal = v;
                best = i;
            }
        }

        if (best == -1)
            throw runtime_error("No valid pivot found");

        order[k] = best;
        used[best] = true;
    }

    double** tmp = new double*[rows];

    for (int i = 0; i < rows; i++) {
        tmp[i] = new double[cols];
        for (int j = 0; j < cols; j++)
            tmp[i][j] = data[order[i]][j];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            data[i][j] = tmp[i][j];
        delete[] tmp[i];
    }
    delete[] tmp;
}