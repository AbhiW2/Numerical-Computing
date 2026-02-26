#include <iostream>
#include "Bisection.hpp"
#include "NewtonRaphson.hpp"
#include "FixedPoint.hpp"

using namespace std;

int main() {

    double a, b, x0;

    // -------- Bisection --------
    cout << "Enter interval [a b] for Bisection: ";
    cin >> a >> b;

    Bisection bis;
    bis.findRoot(a, b);

    // -------- Newton --------
    cout << "\nEnter initial guess for Newton: ";
    cin >> x0;

    RootSolver* solver1 = new NewtonRaphson();
    solver1->findRoot(x0);
    delete solver1;

    // -------- Fixed Point --------
    cout << "\nEnter initial guess for Fixed Point: ";
    cin >> x0;

    RootSolver* solver2 = new FixedPoint();
    solver2->findRoot(x0);
    delete solver2;

    return 0;
}
