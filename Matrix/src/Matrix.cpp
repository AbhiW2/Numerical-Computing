#include "../Include/Matrix.hpp"

using namespace std;

// Constructor
Matrix::Matrix(int r, int c)
{
    if (r <= 0 || c <= 0)
        throw invalid_argument("Matrix dimensions must be positive");

    rows = r;
    cols = c;

    data = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        data[i] = new double[cols];
        for (int j = 0; j < cols; j++)
            data[i][j] = 0.0;
    }
}

// Copy Constructor (Deep Copy)
Matrix::Matrix(const Matrix &other)
{
    rows = other.rows;
    cols = other.cols;

    data = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        data[i] = new double[cols];
        for (int j = 0; j < cols; j++)
            data[i][j] = other.data[i][j];
    }
}

// Destructor
Matrix::~Matrix()
{
    for (int i = 0; i < rows; i++)
        delete[] data[i];

    delete[] data;
}

// Operator +
Matrix Matrix::operator+(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
        throw invalid_argument("Matrix dimensions must match for addition");

    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] + other.data[i][j];

    return result;
}

// Operator -
Matrix Matrix::operator-(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
        throw invalid_argument("Matrix dimensions must match for subtraction");

    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] - other.data[i][j];

    return result;
}

// Operator *
Matrix Matrix::operator*(const Matrix &other) const
{
    if (cols != other.rows)
        throw invalid_argument("Invalid matrix dimensions for multiplication");

    Matrix result(rows, other.cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < other.cols; j++)
            for (int k = 0; k < cols; k++)
                result.data[i][j] += data[i][k] * other.data[k][j];

    return result;
}

// Set value
void Matrix::setValue(int r, int c, double value)
{
    data[r][c] = value;
}

// Get value
double Matrix::getValue(int r, int c) const
{
    return data[r][c];
}

// Display matrix
void Matrix::display() const
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            cout << data[i][j] << " ";
        cout << endl;
    }
}

// Getters
int Matrix::getRows() const
{
    return rows;
}

int Matrix::getCols() const
{
    return cols;
}