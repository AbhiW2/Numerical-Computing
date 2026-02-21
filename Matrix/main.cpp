#include "Matrix.hpp"

int main() {

    int r, c;
    cout << "Enter rows and columns: ";
    cin >> r >> c;

    Matrix m1(r, c), m2(r, c);

    cout << "\nMatrix 1\n";
    m1.input();

    cout << "\nMatrix 2\n";
    m2.input();

    // Addition
    Matrix sum = m1 + m2;
    cout << "\nAddition of two matrices:";
    sum.display();

    // Subtraction
    Matrix diff = m1 - m2;
    cout << "\nSubtraction of two matrices:";
    diff.display();

    // Gaussian Elimination on first matrix
    cout << "\nPerforming Gaussian Elimination on Matrix 1";
    m1.gaussianElimination();

    return 0;
}