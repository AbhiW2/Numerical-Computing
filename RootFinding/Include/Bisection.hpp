

#ifndef BISECTION_HPP
#define BISECTION_HPP

// -------- Parent Class --------
class RootSolver {
public:
    virtual double findRoot(double x0) = 0;  // Pure virtual
    virtual ~RootSolver() {}
};

// -------- Bisection Child Class --------
class Bisection : public RootSolver {
public:
    double func(double x);

    // Required override (not used directly)
    double findRoot(double x0) override;

    double findRoot(double a, double b);
};

#endif
