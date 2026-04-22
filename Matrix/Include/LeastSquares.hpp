#ifndef LEAST_SQUARES_HPP
#define LEAST_SQUARES_HPP

#include "Interpolation.hpp"
#include <vector>

class LeastSquares : public Interpolation {
public:
    LeastSquares(int n, int degree);

    void   loadDataPoints(std::ifstream& fin) override;

    double interpolate(double x_query) override;


    void fit();

    std::vector<double> getCoefficients() const;

    void printFormula() const;

    void printNormalEquations() const;

    std::vector<double> residuals() const;

    // Return the Sum-of-Squared-Errors
    double SSE() const;

private:
    int                 m;       
    std::vector<double> coeffs;  
    bool                fitted;  

    double evalPoly(double xv) const;

    static std::vector<double> solveAugmented(
            std::vector<std::vector<double>> aug, int sz);

    static constexpr double EPS = 1e-14;
};

#endif // LEAST_SQUARES_HPP