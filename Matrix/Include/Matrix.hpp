#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <fstream>
#include <stdexcept>

class Matrix {
protected:
    int rows;
    int cols;
    double** data;

public:
    Matrix();
    Matrix(int rows, int cols);
    Matrix(const Matrix& m);
    virtual ~Matrix();

    int getRows() const { return rows; }
    int getCols() const { return cols; }

    void readFromFile(std::ifstream& fin);
    void displayToFile(std::ofstream& fout) const;

    Matrix operator+(const Matrix& m) const;
    Matrix operator-(const Matrix& m) const;
    Matrix operator*(const Matrix& m) const;

    double&       operator()(int i, int j);
    const double& operator()(int i, int j) const;

    bool operator==(const Matrix& m) const;

    friend std::istream& operator>>(std::istream& in, Matrix& m);
    friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

    bool isSquare()             const;
    bool isSymmetric()          const;
    bool isIdentity()           const;
    bool isNull()               const;
    bool isDiagonal()           const;
    bool isDiagonallyDominant() const;
    bool isTranspose(const Matrix& m) const;

    void   makeDiagonal();
    void   makeDiagonallyDominant();  // greedy row-swap reordering
    Matrix transpose() const;

    double determinant() const;
    Matrix inverse()     const;
};

#endif