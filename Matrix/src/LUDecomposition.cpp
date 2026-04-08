#include "../Include/LUDecomposition.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <numeric>

using namespace std;

LUDecomposition::LUDecomposition(int n, LUMethod m)
    : LinearSystemSolver(n, n + 1), L(n, n), U(n, n), method(m)
{}

void LUDecomposition::setMethod(LUMethod newMethod)
{
    method = newMethod;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < rows; j++) { L(i,j) = 0.0; U(i,j) = 0.0; }
}

void LUDecomposition::loadAugmentedMatrix(ifstream& Afile, ifstream& /*bfile*/)
{
    loadFromAugmentedFile(Afile);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < rows; j++) { L(i,j) = 0.0; U(i,j) = 0.0; }
}

// ── DOOLITTLE ─────────────────────────────────────────────────────────────────
vector<int> LUDecomposition::doolittle()
{
    int n = rows;
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U(i,j) = data[i][j];

    for (int k = 0; k < n; k++) {
        int pivot = k;
        double maxVal = fabs(U(k,k));
        for (int i = k+1; i < n; i++)
            if (fabs(U(i,k)) > maxVal) { maxVal = fabs(U(i,k)); pivot = i; }

        if (pivot != k) {
            for (int j = 0; j < n; j++) swap(U(k,j), U(pivot,j));
            for (int j = 0; j < k; j++) swap(L(k,j), L(pivot,j));
            swap(perm[k], perm[pivot]);
        }

        L(k,k) = 1.0;
        if (fabs(U(k,k)) < EPS) continue;

        for (int i = k+1; i < n; i++) {
            L(i,k) = U(i,k) / U(k,k);
            for (int j = k; j < n; j++)
                U(i,j) -= L(i,k) * U(k,j);
        }
    }
    return perm;
}

// ── CROUT ─────────────────────────────────────────────────────────────────────
vector<int> LUDecomposition::crout()
{
    int n = rows;
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);

    Matrix A(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A(i,j) = data[i][j];

    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) sum += L(i,k) * U(k,j);
            L(i,j) = A(i,j) - sum;
        }

        int pivot = j;
        double maxVal = fabs(L(j,j));
        for (int i = j+1; i < n; i++)
            if (fabs(L(i,j)) > maxVal) { maxVal = fabs(L(i,j)); pivot = i; }

        if (pivot != j) {
            for (int k = 0; k < n; k++) swap(L(j,k), L(pivot,k));
            for (int k = 0; k < n; k++) swap(A(j,k), A(pivot,k));
            for (int k = 0; k < j; k++) swap(U(j,k), U(pivot,k));
            swap(perm[j], perm[pivot]);
        }

        U(j,j) = 1.0;
        if (fabs(L(j,j)) < EPS) continue;

        for (int i = j+1; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) sum += L(j,k) * U(k,i);
            U(j,i) = (A(j,i) - sum) / L(j,j);
        }
    }
    return perm;
}

// ── CHOLESKY ──────────────────────────────────────────────────────────────────
vector<int> LUDecomposition::cholesky()
{
    int n = rows;
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (fabs(data[i][j] - data[j][i]) > EPS)
                throw runtime_error("Cholesky requires a symmetric matrix");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) sum += L(i,k) * L(j,k);
            if (i == j) {
                double val = data[i][i] - sum;
                if (val <= EPS)
                    throw runtime_error(
                        "Cholesky: not positive-definite at diagonal " + to_string(i));
                L(i,j) = sqrt(val);
            } else {
                if (fabs(L(j,j)) < EPS)
                    throw runtime_error(
                        "Cholesky: zero diagonal in L at column " + to_string(j));
                L(i,j) = (data[i][j] - sum) / L(j,j);
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U(i,j) = L(j,i);

    return perm;
}

vector<int> LUDecomposition::decompose()
{
    switch (method) {
        case LUMethod::DOOLITTLE: return doolittle();
        case LUMethod::CROUT:     return crout();
        case LUMethod::CHOLESKY:  return cholesky();
        default: throw runtime_error("Unknown LU method");
    }
}

vector<double> LUDecomposition::forwardSubstitution(const Matrix& Lmat,
                                                     const vector<double>& b)
{
    int n = rows;
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++) {
        if (fabs(Lmat(i,i)) < EPS) continue;
        double sum = 0.0;
        for (int j = 0; j < i; j++) sum += Lmat(i,j) * y[j];
        y[i] = (b[i] - sum) / Lmat(i,i);
    }
    return y;
}

vector<double> LUDecomposition::backSubstitution(const Matrix& Umat,
                                                  const vector<double>& y)
{
    int n = rows;
    vector<double> x(n, 0.0);
    for (int i = n-1; i >= 0; i--) {
        if (fabs(Umat(i,i)) < EPS) continue;
        double sum = 0.0;
        for (int j = i+1; j < n; j++) sum += Umat(i,j) * x[j];
        x[i] = (y[i] - sum) / Umat(i,i);
    }
    return x;
}

vector<double> LUDecomposition::solve()
{
    int n = rows;

    vector<double> b_orig(n);
    for (int i = 0; i < n; i++) b_orig[i] = data[i][cols - 1];

    vector<int> perm = decompose();

    vector<double> b_perm(n);
    for (int i = 0; i < n; i++) b_perm[i] = b_orig[perm[i]];

    int rank = 0;
    for (int i = 0; i < n; i++) if (fabs(U(i,i)) > EPS) rank++;

    if (rank < n)
        cout << "Note: rank-deficient — returning particular solution (free variables = 0).\n\n";
    else
        cout << "System has a unique solution.\n\n";

    vector<double> y = forwardSubstitution(L, b_perm);
    return backSubstitution(U, y);
}