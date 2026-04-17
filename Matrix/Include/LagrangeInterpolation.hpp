#ifndef LAGRANGE_INTERPOLATION_HPP
#define LAGRANGE_INTERPOLATION_HPP

#include "Interpolation.hpp"
#include <vector>

class LagrangeInterpolation : public Interpolation {
public:
    explicit LagrangeInterpolation(int n);

    void loadDataPoints(std::ifstream& fin) override;

    double interpolate(double x_query) override;

    std::vector<double> basisPolynomials(double x_query) const;

    void printFormula() const;


private:
    double basisAt(int i, double x) const;

    static constexpr double EPS = 1e-14;
};

#endif 