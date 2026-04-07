#include "../Include/Matrix.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <vector>

using namespace std;

Matrix::Matrix() : rows(0), cols(0), data(nullptr) {}

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(nullptr)
{
    if (r <= 0 || c <= 0)
        throw invalid_argument("Matrix dimensions must be positive");
    data = new double*[rows];
    for (int i = 0; i < rows; ++i)
        data[i] = new double[cols]();
}

Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(nullptr)
{
    data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; ++j)
            data[i][j] = other.data[i][j];
    }
}

Matrix::~Matrix()
{
    if (data) {
        for (int i = 0; i < rows; ++i) delete[] data[i];
        delete[] data;
        data = nullptr;
    }
}

void Matrix::readFromFile(ifstream& fin)
{
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            fin >> data[i][j];
}

void Matrix::displayToFile(ofstream& fout) const { fout << *this; }

Matrix Matrix::operator+(const Matrix& m) const
{
    if (rows != m.rows || cols != m.cols)
        throw invalid_argument("Dimension mismatch for addition");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result(i,j) = (*this)(i,j) + m(i,j);
    return result;
}

Matrix Matrix::operator-(const Matrix& m) const
{
    if (rows != m.rows || cols != m.cols)
        throw invalid_argument("Dimension mismatch for subtraction");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result(i,j) = (*this)(i,j) - m(i,j);
    return result;
}

Matrix Matrix::operator*(const Matrix& m) const
{
    if (cols != m.rows)
        throw invalid_argument("Invalid dimensions for multiplication");
    Matrix result(rows, m.cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < m.cols; ++j)
            for (int k = 0; k < cols; ++k)
                result(i,j) += (*this)(i,k) * m(k,j);
    return result;
}

double& Matrix::operator()(int i, int j)
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw out_of_range("Index out of range");
    return data[i][j];
}

const double& Matrix::operator()(int i, int j) const
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw out_of_range("Index out of range");
    return data[i][j];
}

bool Matrix::operator==(const Matrix& m) const
{
    if (rows != m.rows || cols != m.cols) return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (fabs((*this)(i,j) - m(i,j)) > 1e-10) return false;
    return true;
}

istream& operator>>(istream& in, Matrix& m)
{
    for (int i = 0; i < m.rows; ++i)
        for (int j = 0; j < m.cols; ++j)
            in >> m(i,j);
    return in;
}

ostream& operator<<(ostream& out, const Matrix& m)
{
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j)
            out << setw(12) << fixed << setprecision(6) << m(i,j) << " ";
        out << "\n";
    }
    return out;
}

bool Matrix::isSquare()     const { return rows == cols; }

bool Matrix::isSymmetric() const
{
    if (!isSquare()) return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (fabs((*this)(i,j) - (*this)(j,i)) > 1e-10) return false;
    return true;
}

bool Matrix::isIdentity() const
{
    if (!isSquare()) return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs((*this)(i,j) - expected) > 1e-10) return false;
        }
    return true;
}

bool Matrix::isNull() const
{
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (fabs((*this)(i,j)) > 1e-10) return false;
    return true;
}

bool Matrix::isDiagonal() const
{
    if (!isSquare()) return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (i != j && fabs((*this)(i,j)) > 1e-10) return false;
    return true;
}

bool Matrix::isDiagonallyDominant() const
{
    if (!isSquare()) return false;
    for (int i = 0; i < rows; ++i) {
        double off_sum = 0.0;
        for (int j = 0; j < cols; ++j)
            if (i != j) off_sum += fabs((*this)(i,j));
        if (fabs((*this)(i,i)) <= off_sum) return false;
    }
    return true;
}

bool Matrix::isTranspose(const Matrix& m) const
{
    if (rows != m.cols || cols != m.rows) return false;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (fabs((*this)(i,j) - m(j,i)) > 1e-10) return false;
    return true;
}

void Matrix::makeDiagonal()
{
    if (!isSquare())
        throw invalid_argument("makeDiagonal requires square matrix");
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            (*this)(i,j) = (i == j) ? 1.0 : 0.0;
}

// Called ONLY by GaussJacobi and GaussSeidel when isDiagonallyDominant() == false.
// Greedy column-maximum row reordering: for each column k, swap in the
// unused row that has the largest |value| in that column.
// If strict dominance still cannot be achieved (rank-deficient matrix),
// a warning is printed and the best reordering is kept.
void Matrix::makeDiagonallyDominant()
{
    int n = rows;   // square part size (cols may be n+1 for augmented)

    vector<bool> used(n, false);
    vector<int>  perm(n, -1);

    for (int col = 0; col < n; ++col) {
        int    best     = -1;
        double best_val = -1.0;
        for (int row = 0; row < n; ++row) {
            if (!used[row] && fabs(data[row][col]) > best_val) {
                best_val = fabs(data[row][col]);
                best     = row;
            }
        }
        perm[col]    = best;
        used[best]   = true;
    }

    // Apply permutation via temp buffer
    double** tmp = new double*[n];
    for (int i = 0; i < n; ++i) {
        tmp[i] = new double[cols];
        for (int j = 0; j < cols; ++j)
            tmp[i][j] = data[perm[i]][j];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < cols; ++j)
            data[i][j] = tmp[i][j];
        delete[] tmp[i];
    }
    delete[] tmp;

    if (isDiagonallyDominant()) {
        cout << "Matrix is now diagonally dominant after row reordering.\n";
    } else {
        cout << "Warning: strict diagonal dominance could not be achieved\n"
             << "         (matrix may be rank-deficient).\n"
             << "         Best achievable reordering has been applied.\n"
             << "         Gaussian elimination / LU still work correctly.\n"
             << "         Iterative methods (Jacobi/Seidel) may not converge.\n";
    }
}

Matrix Matrix::transpose() const
{
    Matrix t(cols, rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            t(j,i) = (*this)(i,j);
    return t;
}

double Matrix::determinant() const
{
    if (rows != cols) throw runtime_error("Determinant only for square matrix");
    int n = rows;
    Matrix temp = *this;
    double det  = 1.0;
    for (int k = 0; k < n; k++) {
        int pivot = k;
        for (int i = k+1; i < n; i++)
            if (fabs(temp(i,k)) > fabs(temp(pivot,k))) pivot = i;
        if (fabs(temp(pivot,k)) < 1e-12) return 0.0;
        if (pivot != k) {
            for (int j = 0; j < n; j++) swap(temp(k,j), temp(pivot,j));
            det *= -1;
        }
        det *= temp(k,k);
        for (int i = k+1; i < n; i++) {
            double factor = temp(i,k) / temp(k,k);
            for (int j = k; j < n; j++)
                temp(i,j) -= factor * temp(k,j);
        }
    }
    return det;
}

Matrix Matrix::inverse() const
{
    if (rows != cols) throw runtime_error("Inverse only for square matrix");
    int n = rows;
    Matrix aug(n, 2*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug(i,j) = (*this)(i,j);
        aug(i, n+i) = 1.0;
    }
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int k = i+1; k < n; k++)
            if (fabs(aug(k,i)) > fabs(aug(pivot,i))) pivot = k;
        if (fabs(aug(pivot,i)) < 1e-12)
            throw runtime_error("Matrix not invertible");
        if (pivot != i)
            for (int j = 0; j < 2*n; j++) swap(aug(i,j), aug(pivot,j));
        double d = aug(i,i);
        for (int j = 0; j < 2*n; j++) aug(i,j) /= d;
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = aug(k,i);
            for (int j = 0; j < 2*n; j++) aug(k,j) -= factor * aug(i,j);
        }
    }
    Matrix inv(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv(i,j) = aug(i, j+n);
    return inv;
}