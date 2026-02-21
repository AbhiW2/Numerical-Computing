

#include <iostream>
#include <cmath>
#include "NewtonRaphson.hpp"

using namespace std;

double NewtonRaphson::func(double x) {
    return 4 * x * x * x - 3 * x;
}

double NewtonRaphson::derivative(double x) {
    return 12 * x * x - 3;
}

double NewtonRaphson::findRoot(double x0) {

    double tolerance = 0.0001;
    int iterations = 0;
    double x1;

    do {
        if (derivative(x0) == 0) {
            cout << "Derivative zero. Method fails.\n";
            return -1;
        }

        x1 = x0 - func(x0) / derivative(x0);
        x0 = x1;
        iterations++;

    } while (abs(func(x1)) > tolerance);

    cout << "Newton Root: " << x1 << endl;
    cout << "Iterations: " << iterations << endl;

    return x1;
}
