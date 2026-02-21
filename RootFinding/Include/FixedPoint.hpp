
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP

#include "Bisection.hpp"   // IMPORTANT

class FixedPoint : public RootSolver {
public:
    double g(double x);
    double findRoot(double x0) override;
};

#endif
