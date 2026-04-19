#include "../Include/LagrangeInterpolation.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

using namespace std;


LagrangeInterpolation::LagrangeInterpolation(int n)
    : Interpolation(n)
{}


void LagrangeInterpolation::loadDataPoints(ifstream& fin) {
    Interpolation::loadDataPoints(fin);

    // Guard: all x nodes must be distinct
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (fabs(x(i) - x(j)) < EPS)
                throw runtime_error(
                    "Lagrange requires distinct x nodes: "
                    "x[" + to_string(i) + "] == x[" + to_string(j) + "]");
}

double LagrangeInterpolation::basisAt(int i, double xq) const {
    double result = 1.0;
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        double denom = x(i) - x(j);
        if (fabs(denom) < EPS)
            throw runtime_error("Duplicate x nodes detected during evaluation.");
        result *= (xq - x(j)) / denom;
    }
    return result;
}

double LagrangeInterpolation::interpolate(double x_query) {
    double P = 0.0;
    for (int i = 0; i < n; ++i)
        P += y(i) * basisAt(i, x_query);
    return P;
}


vector<double> LagrangeInterpolation::basisPolynomials(double x_query) const {
    vector<double> L(n);
    for (int i = 0; i < n; ++i)
        L[i] = basisAt(i, x_query);
    return L;
}

void LagrangeInterpolation::printFormula() const {
    cout << fixed << setprecision(4);
    cout << "Lagrange Polynomial P(x):\n";
    cout << "=========================\n\n";

    for (int i = 0; i < n; ++i) {
        cout << "  Term " << i << ":  " << y(i) << " * L_" << i << "(x)\n";
        cout << "         where  L_" << i << "(x) = ";

        bool firstFactor = true;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double denom = x(i) - x(j);
            if (!firstFactor) cout << " * ";

            // print  (x - x_j) / (x_i - x_j)
            cout << "(x";
            if (x(j) >= 0) cout << " - " <<  x(j);
            else            cout << " + " << -x(j);
            cout << ") / " << denom;

            firstFactor = false;
        }
        cout << "\n\n";
    }
}

