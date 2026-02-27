// #ifndef MATRIX_HPP
// #define MATRIX_HPP

// #include <iostream>
// using namespace std;

// class Matrix {
// protected:
//     int rows, cols;
//     double **data;

// public:
//     Matrix(int r, int c);
//     Matrix(const Matrix &other);
//     virtual ~Matrix();

//     // Basic operations
//     Matrix add(const Matrix &other);
//     Matrix subtract(const Matrix &other);
//     Matrix multiply(const Matrix &other);

//     void setValue(int r, int c, double val);
//     double getValue(int r, int c) const;

//     void display() const;
// };

// #endif




#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>

class Matrix {
protected:
    int rows, cols;
    double **data;

public:
    Matrix(int r, int c);
    Matrix(const Matrix &other);
    virtual ~Matrix();

    Matrix add(const Matrix &other);
    Matrix subtract(const Matrix &other);
    Matrix multiply(const Matrix &other);

    void setValue(int r, int c, double val);
    double getValue(int r, int c) const;

    void display() const;
};

#endif