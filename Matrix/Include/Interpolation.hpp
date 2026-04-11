#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "Matrix.hpp"
#include <vector>
#include <fstream>

class Interpolation : public Matrix {
public:
    Interpolation(int n);           // n = number of data points  →  Matrix(n, 2)
    virtual ~Interpolation() = default;

    virtual void loadDataPoints(std::ifstream& fin);

    virtual double interpolate(double x_query) = 0;
    
    void displayDataPoints() const;

protected:
    int n;   

    double  x(int i) const { return data[i][0]; }
    double  y(int i) const { return data[i][1]; }
    double& x(int i)       { return data[i][0]; }
    double& y(int i)       { return data[i][1]; }
};

#endif 