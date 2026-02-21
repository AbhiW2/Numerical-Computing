#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
using namespace std;


class MatrixBase {
public:
    virtual void input() = 0;
    virtual void display() const = 0;
    virtual void gaussianElimination() = 0;
};

/*
   Derived Class (Inheritance)
   Encapsulation: data members are private
  */
class Matrix : public MatrixBase {
private:
    int rows, cols;
    double data[10][10];   // fixed size for simplicity

public:
    // Constructor
    Matrix(int r = 0, int c = 0);

    // Operator Overloading
    Matrix operator+(const Matrix& other);
    Matrix operator-(const Matrix& other);

    // Overriding Base Class Functions (Polymorphism)
    void input() override;
    void display() const override;
    void gaussianElimination() override;
};

#endif