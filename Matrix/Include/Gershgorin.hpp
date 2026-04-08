#ifndef GERSHGORIN_HPP
#define GERSHGORIN_HPP

#include "Matrix.hpp"
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

class Gershgorin : public Matrix {
private:
    int   n;
    std::vector<double>  centers;
    std::vector<double>  radii;
    std::vector<double>  lowers;
    std::vector<double>  uppers;

    void compute();

public:
    explicit Gershgorin(const Matrix& matrix);

    // Core 
    std::pair<double, double> eigenvalueBounds()    const;
    double                    spectralRadiusBound() const;

    // Properties
    bool containsZero()         const;
    bool allDiscsPositive()     const;
    bool isDiagonallyDominant() const;

    // Display
    void print()        const;
    void printSummary() const;
};

#endif