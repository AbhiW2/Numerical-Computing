
#include "../Include/GaussianElimination.hpp"
#include <cmath>
#include <iostream>

using namespace std;

static const double EPS = 1e-12;

// Constructor
GaussianElimination::GaussianElimination(int r, int c)
    : Matrix(r, c) {}

void GaussianElimination::loadAugmentedMatrix(ifstream &Afile, ifstream &bfile)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols - 1; j++)
            Afile >> data[i][j];

    for (int i = 0; i < rows; i++)
        bfile >> data[i][cols - 1];
}

// Forward elimination with pivot normalization
void GaussianElimination::forwardElimination()
{
    for (int i = 0; i < rows; i++)
    {
        // partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < rows; k++)
            if (fabs(data[k][i]) > fabs(data[maxRow][i]))
                maxRow = k;

        if (fabs(data[maxRow][i]) < EPS)
            continue;

        if (maxRow != i)
            for (int j = 0; j < cols; j++)
                swap(data[i][j], data[maxRow][j]);

        // Normalize pivot row
        double pivot = data[i][i];
        for (int j = i; j < cols; j++)
            data[i][j] /= pivot;

        // Eliminate below
        for (int k = i + 1; k < rows; k++)
        {
            double factor = data[k][i];
            for (int j = i; j < cols; j++)
                data[k][j] -= factor * data[i][j];
        }
    }
}

// Back substitution (no division needed)
vector<double> GaussianElimination::backSubstitution()
{
    vector<double> x(rows, 0.0);

    for (int i = rows - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < rows; j++)
            sum += data[i][j] * x[j];

        x[i] = data[i][cols - 1] - sum;
    }

    return x;
}

// Solve Ax = b
vector<double> GaussianElimination::solve()
{
    forwardElimination();
    return backSubstitution();
}