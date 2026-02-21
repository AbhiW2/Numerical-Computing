#include "Matrix.hpp"

/* Constructor */
Matrix::Matrix(int r, int c) {
    rows = r;
    cols = c;

    // Initialize matrix with 0
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            data[i][j] = 0;
}

/* Input Function */
void Matrix::input() {
    cout << "Enter matrix elements:\n";
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            cin >> data[i][j];
}

/* Display Function */
void Matrix::display() const {
    cout << "\nMatrix:\n";
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            cout << data[i][j] << " ";
        cout << endl;
    }
}

/* Operator Overloading for Addition */
Matrix Matrix::operator+(const Matrix& other) {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] + other.data[i][j];

    return result;
}

/* Operator Overloading for Subtraction */
Matrix Matrix::operator-(const Matrix& other) {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] - other.data[i][j];

    return result;
}

/*
   Gaussian Elimination 
   Converts matrix into upper triangular form
*/
void Matrix::gaussianElimination() {
    if (rows != cols) {
        cout << "Gaussian Elimination requires square matrix.\n";
        return;
    }

    for (int k = 0; k < rows - 1; k++) {
        for (int i = k + 1; i < rows; i++) {

            if (data[k][k] == 0) {
                cout << "Zero pivot encountered. Cannot proceed.\n";
                return;
            }

            double factor = data[i][k] / data[k][k];

            for (int j = k; j < cols; j++) {
                data[i][j] -= factor * data[k][j];
            }
        }
    }

    cout << "\nUpper Triangular Matrix:\n";
    display();
}