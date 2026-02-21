
#ifndef NEWTONRAPHSON_HPP
#define NEWTONRAPHSON_HPP

#include "Bisection.hpp"   // IMPORTANT

class NewtonRaphson : public RootSolver {
public:
    double func(double x);
    double derivative(double x);
    double findRoot(double x0) override;
};

#endif

