#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <stdexcept>

class Matrix
{
protected:
    int rows;
    int cols;
    double **data;

public:
    // Constructor
    Matrix(int r, int c);

    // Copy Constructor (Deep Copy)
    Matrix(const Matrix &other);

    // Destructor
    virtual ~Matrix();

    // Operator Overloading
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;
    Matrix operator*(const Matrix &other) const;

    // Utility Functions
    void setValue(int r, int c, double value);
    double getValue(int r, int c) const;
    void display() const;

    int getRows() const;
    int getCols() const;
};

#endif