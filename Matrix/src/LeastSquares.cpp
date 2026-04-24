#include "../Include/LeastSquares.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <string>

using namespace std;

LeastSquares::LeastSquares(int n, int degree)
    : Interpolation(n),     // calls Matrix(n, 2) through the chain
      m(degree),
      coeffs(degree + 1, 0.0),
      fitted(false)
{
    if (degree < 1)
        throw invalid_argument("LeastSquares: polynomial degree must be >= 1.");
    if (degree >= n)
        throw invalid_argument(
            "LeastSquares: degree must be strictly less than n (number of points).");
}

// ─────────────────────────────────────────────────────────────────────────────
//  loadDataPoints  –  delegates to parent then triggers fit()
// ─────────────────────────────────────────────────────────────────────────────

void LeastSquares::loadDataPoints(ifstream& fin) {
    // Reuse parent's file-reading logic (validates count, fills data[][])
    Interpolation::loadDataPoints(fin);

    // Immediately compute the best-fit coefficients
    fit();
}

// ─────────────────────────────────────────────────────────────────────────────
//  fit  –  build and solve the normal equations  (AᵀA) c = Aᵀb
// ─────────────────────────────────────────────────────────────────────────────

void LeastSquares::fit() {
    int sz = m + 1;   // number of unknowns / size of normal-equation system

    // ── build the (m+1)×(m+1) augmented matrix ──────────────────────────
    // AᵀA  is symmetric;  element (r,c) = Σ x_i^(r+c)
    // Aᵀb  element (r)    = Σ x_i^r * y_i
    vector<vector<double>> aug(sz, vector<double>(sz + 1, 0.0));

    for (int r = 0; r < sz; ++r) {
        for (int c = 0; c < sz; ++c) {
            double sum = 0.0;
            for (int i = 0; i < n; ++i)
                sum += pow(x(i), r + c);
            aug[r][c] = sum;
        }
        // right-hand side
        double rhs = 0.0;
        for (int i = 0; i < n; ++i)
            rhs += pow(x(i), r) * y(i);
        aug[r][sz] = rhs;
    }

    // ── solve the system ─────────────────────────────────────────────────
    coeffs = solveAugmented(aug, sz);
    fitted = true;

    cout << "Least-Squares fit complete (degree " << m
         << ", n=" << n << " points).\n\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  interpolate  –  evaluate P(x_query)          [overrides pure virtual]
// ─────────────────────────────────────────────────────────────────────────────

double LeastSquares::interpolate(double x_query) {
    if (!fitted)
        throw runtime_error("LeastSquares: call loadDataPoints() or fit() first.");
    return evalPoly(x_query);
}

// ─────────────────────────────────────────────────────────────────────────────
//  getCoefficients
// ─────────────────────────────────────────────────────────────────────────────

vector<double> LeastSquares::getCoefficients() const {
    if (!fitted)
        throw runtime_error("LeastSquares: coefficients not yet computed.");
    return coeffs;
}

// ─────────────────────────────────────────────────────────────────────────────
//  printFormula  –  mirrors LagrangeInterpolation::printFormula()
// ─────────────────────────────────────────────────────────────────────────────

void LeastSquares::printFormula() const {
    if (!fitted)
        throw runtime_error("LeastSquares: no formula yet – call fit() first.");

    cout << fixed << setprecision(6);
    cout << "Least-Squares Polynomial P(x) [degree " << m << "]:\n";
    cout << "=============================================\n\n";
    cout << "  P(x) = ";

    bool firstTerm = true;
    for (int k = 0; k <= m; ++k) {
        double c = coeffs[k];
        if (fabs(c) < EPS) continue;

        if (!firstTerm) {
            cout << (c >= 0 ? "  +  " : "  -  ");
            c = fabs(c);
        }

        cout << c;

        if (k == 1)      cout << " * x";
        else if (k > 1)  cout << " * x^" << k;

        firstTerm = false;
    }
    if (firstTerm) cout << "0";
    cout << "\n\n";

    cout << "  Coefficients:\n";
    for (int k = 0; k <= m; ++k)
        cout << "    c[" << k << "] = " << setw(14) << coeffs[k] << "\n";
    cout << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  printNormalEquations  –  educational display
// ─────────────────────────────────────────────────────────────────────────────

void LeastSquares::printNormalEquations() const {
    if (n == 0)
        throw runtime_error("LeastSquares: load data points first.");

    int sz = m + 1;
    cout << fixed << setprecision(4);
    cout << "Normal Equations  (AᵀA) c = Aᵀb:\n";
    cout << "===================================\n";

    vector<vector<double>> aug(sz, vector<double>(sz + 1, 0.0));
    for (int r = 0; r < sz; ++r) {
        for (int c = 0; c < sz; ++c) {
            double sum = 0.0;
            for (int i = 0; i < n; ++i) sum += pow(x(i), r + c);
            aug[r][c] = sum;
        }
        double rhs = 0.0;
        for (int i = 0; i < n; ++i) rhs += pow(x(i), r) * y(i);
        aug[r][sz] = rhs;
    }

    for (int r = 0; r < sz; ++r) {
        cout << "  [ ";
        for (int c = 0; c <= sz; ++c) {
            if (c == sz) cout << " | ";
            cout << setw(12) << aug[r][c];
        }
        cout << " ]\n";
    }
    cout << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  residuals
// ─────────────────────────────────────────────────────────────────────────────

vector<double> LeastSquares::residuals() const {
    if (!fitted)
        throw runtime_error("LeastSquares: call fit() first.");

    vector<double> res(n);
    for (int i = 0; i < n; ++i)
        res[i] = y(i) - evalPoly(x(i));
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
//  SSE
// ─────────────────────────────────────────────────────────────────────────────

double LeastSquares::SSE() const {
    auto r = residuals();
    double s = 0.0;
    for (double ri : r) s += ri * ri;
    return s;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Private helpers
// ─────────────────────────────────────────────────────────────────────────────

double LeastSquares::evalPoly(double xv) const {
    // Horner's method for numerical stability
    double result = coeffs[m];
    for (int k = m - 1; k >= 0; --k)
        result = result * xv + coeffs[k];
    return result;
}

vector<double> LeastSquares::solveAugmented(
        vector<vector<double>> aug, int sz)
{
    // Forward elimination with partial pivoting
    for (int col = 0; col < sz; ++col) {
        int pivot = col;
        for (int row = col + 1; row < sz; ++row)
            if (fabs(aug[row][col]) > fabs(aug[pivot][col]))
                pivot = row;
        swap(aug[col], aug[pivot]);

        double diag = aug[col][col];
        if (fabs(diag) < EPS)
            throw runtime_error(
                "LeastSquares: singular normal-equation matrix "
                "(try a lower polynomial degree).");

        for (int row = col + 1; row < sz; ++row) {
            double factor = aug[row][col] / diag;
            for (int k = col; k <= sz; ++k)
                aug[row][k] -= factor * aug[col][k];
        }
    }

    // Back-substitution
    vector<double> sol(sz, 0.0);
    for (int row = sz - 1; row >= 0; --row) {
        sol[row] = aug[row][sz];
        for (int k = row + 1; k < sz; ++k)
            sol[row] -= aug[row][k] * sol[k];
        sol[row] /= aug[row][row];
    }
    return sol;
}