#include "../Include/Gershgorin.hpp"

using namespace std;

Gershgorin::Gershgorin(const Matrix& matrix)
    : Matrix(matrix), n(matrix.getRows())
{
    if (!isSquare())
        throw runtime_error("Gershgorin: matrix must be square");
    compute();
}

void Gershgorin::compute()
{
    centers.resize(n);
    radii  .resize(n);
    lowers .resize(n);
    uppers .resize(n);

    for (int i = 0; i < n; i++) {
        centers[i] = data[i][i];
        radii[i]   = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) radii[i] += fabs(data[i][j]);

        lowers[i] = centers[i] - radii[i];
        uppers[i] = centers[i] + radii[i];
    }
}

pair<double, double> Gershgorin::eigenvalueBounds() const
{
    double lo =  1e300;
    double hi = -1e300;
    for (int i = 0; i < n; i++) {
        if (lowers[i] < lo) lo = lowers[i];
        if (uppers[i] > hi) hi = uppers[i];
    }
    return { lo, hi };
}

double Gershgorin::spectralRadiusBound() const
{
    double maxVal = 0.0;
    for (int i = 0; i < n; i++) {
        double v = fabs(centers[i]) + radii[i];
        if (v > maxVal) maxVal = v;
    }
    return maxVal;
}

bool Gershgorin::containsZero() const
{
    for (int i = 0; i < n; i++)
        if (fabs(centers[i]) <= radii[i]) return true;
    return false;
}

bool Gershgorin::allDiscsPositive() const
{
    for (int i = 0; i < n; i++)
        if (lowers[i] <= 0.0) return false;
    return true;
}

bool Gershgorin::isDiagonallyDominant() const
{
    for (int i = 0; i < n; i++)
        if (fabs(centers[i]) <= radii[i]) return false;
    return true;
}

void Gershgorin::print() const
{
    cout << "\n========== Gershgorin's Theorem ==========\n";
    cout << setw(5)  << "Row"
         << setw(14) << "Center"
         << setw(14) << "Radius"
         << setw(14) << "Lower"
         << setw(14) << "Upper"
         << setw(10) << "Has 0?" << "\n";
    cout << string(71, '-') << "\n";

    for (int i = 0; i < n; i++) {
        bool hasZero = (fabs(centers[i]) <= radii[i]);
        cout << setw(5)  << i + 1
             << setw(14) << fixed << setprecision(4) << centers[i]
             << setw(14) << radii[i]
             << setw(14) << lowers[i]
             << setw(14) << uppers[i]
             << setw(10) << (hasZero ? "YES" : "no") << "\n";
    }

    printSummary();
}

void Gershgorin::printSummary() const
{
    auto [lo, hi] = eigenvalueBounds();

    cout << "\n---------- Summary ----------\n";
    cout << "Eigenvalue bound interval : [ "
         << fixed << setprecision(4) << lo << " , " << hi << " ]\n";
    cout << "Spectral radius bound     : "
         << spectralRadiusBound() << "\n";
    cout << "Any disc contains 0       : "
         << (containsZero()         ? "YES - matrix may be singular"         : "No") << "\n";
    cout << "All discs strictly > 0    : "
         << (allDiscsPositive()     ? "YES - likely positive definite"        : "No") << "\n";
    cout << "Diagonally dominant       : "
         << (isDiagonallyDominant() ? "YES - iterative methods will converge" : "No") << "\n";
}