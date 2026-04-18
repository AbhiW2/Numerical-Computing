#include "../Include/Interpolation.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

using namespace std;

Interpolation::Interpolation(int n)
    : Matrix(n, 2), n(n)
{}

void Interpolation::loadDataPoints(ifstream& fin) {
    int file_n;
    if (!(fin >> file_n))
        throw runtime_error("Cannot read number of data points from file.");

    if (file_n != n)
        throw runtime_error(
            "Size mismatch: Interpolation expects " + to_string(n) +
            " data points but file has " + to_string(file_n) + ".");

    for (int i = 0; i < n; ++i) {
        if (!(fin >> data[i][0] >> data[i][1]))
            throw runtime_error(
                "Error reading data point " + to_string(i) + " from file.");
    }

    cout << "Data points loaded (" << n << " points).\n\n";
}

void Interpolation::displayDataPoints() const {
    cout << fixed << setprecision(6);
    cout << "  i       x_i            y_i\n";
    cout << "-----  ----------     ----------\n";
    for (int i = 0; i < n; ++i)
        cout << "  " << i << "    " << setw(12) << data[i][0]
             << "    " << setw(12) << data[i][1] << "\n";
    cout << "\n";
}