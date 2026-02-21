
#include <iostream>
#include <cmath>
#include "FixedPoint.hpp"

using namespace std;

double FixedPoint::g(double x) {
    return sqrt(3.0 / (4.0 * x));
}

double FixedPoint::findRoot(double x0) {

    double tolerance = 0.0001;
    int iterations = 0;
    double x1;

    do {
        x1 = g(x0);
        x0 = x1;
        iterations++;

    } while (abs(g(x1) - x1) > tolerance);

    cout << "Fixed Point Root: " << x1 << endl;
    cout << "Iterations: " << iterations << endl;

    return x1;
}
