
#include <iostream>
#include <cmath>
#include "Bisection.hpp"

using namespace std;

double Bisection::func(double x) {
    return 4 * x * x * x - 3 * x;
}

// Dummy override (required because parent has pure virtual)
double Bisection::findRoot(double x0) {
    return 0;
}

double Bisection::findRoot(double a, double b) {

    double tolerance = 0.0001;
    int iterations = 0;
    double c;

    if (func(a) * func(b) >= 0) {
        cout << "Invalid interval!\n";
        return -1;
    }

    while ((b - a) >= tolerance) {
        c = (a + b) / 2;

        if (func(c) == 0.0)
            break;
        else if (func(a) * func(c) < 0)
            b = c;
        else
            a = c;

        iterations++;
    }

    cout << "Bisection Root: " << c << endl;
    cout << "Iterations: " << iterations << endl;

    return c;
}
