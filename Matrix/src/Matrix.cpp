
#include "../Include/Matrix.hpp"
#include <iostream>

// Constructor
Matrix::Matrix(int r, int c) {
    rows = r;
    cols = c;

    data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            data[i][j] = 0.0;
        }
    }
}

// Copy Constructor
Matrix::Matrix(const Matrix &other) {
    rows = other.rows;
    cols = other.cols;

    data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            data[i][j] = other.data[i][j];
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for (int i = 0; i < rows; i++) {
        delete[] data[i];
    }
    delete[] data;
}

// Addition
Matrix Matrix::add(const Matrix &other) {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }

    return result;
}

// Subtraction
Matrix Matrix::subtract(const Matrix &other) {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }

    return result;
}

// Multiplication
Matrix Matrix::multiply(const Matrix &other) {
    Matrix result(rows, other.cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < other.cols; j++) {
            for (int k = 0; k < cols; k++) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }

    return result;
}

// Set value
void Matrix::setValue(int r, int c, double val) {
    data[r][c] = val;
}

// Get value
double Matrix::getValue(int r, int c) const {
    return data[r][c];
}

// Display matrix
void Matrix::display() const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}