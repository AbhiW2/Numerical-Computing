#include <iostream>
#include <fstream>
#include "Include/GaussianElimination.hpp"
#include "Include/LUDecomposition.hpp"
#include "Include/GaussJacobi.hpp"
#include "Include/GaussSeidel.hpp"
#include "Include/Matrix.hpp"
#include "Include/Gershgorin.hpp"
#include "Include/LagrangeInterpolation.hpp"
#include "Include/LeastSquares.hpp"

using namespace std;

int main() {
    int mainChoice;

    cout << "===== MENU =====\n";
    cout << "1. Solve Linear System\n";
    cout << "2. Matrix Operations\n";
    cout << "3. Interpolation\n";
    cout << "Enter choice: ";
    cin >> mainChoice;

    if (mainChoice == 1) {
        ifstream peek("49l.txt");
        if (!peek) { cout << "Error: cannot open 49l.txt\n"; return 1; }
        int n, m;
        peek >> n >> m;
        peek.close();
        cout << "\nDetected " << n << " equations from 49l.txt.\n";

        cout << "\nChoose Method:\n";
        cout << "1. Gaussian Elimination\n";
        cout << "2. LU Decomposition (Doolittle)\n";
        cout << "3. LU Decomposition (Crout)\n";
        cout << "4. LU Decomposition (Cholesky)\n";
        cout << "5. Gauss-Jacobi\n";
        cout << "6. Gauss-Seidel\n";
        cout << "Enter choice: ";

        int choice;
        cin >> choice;

        ifstream Afile("49l.txt");
        ifstream bfile("49r.txt");

        if (!Afile) { cout << "Error opening 49l.txt\n"; return 1; }

        try {
            vector<double> result;

            if (choice == 1) {
                GaussianElimination solver(n);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else if (choice == 2) {
                LUDecomposition solver(n, LUMethod::DOOLITTLE);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else if (choice == 3) {
                LUDecomposition solver(n, LUMethod::CROUT);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else if (choice == 4) {
                LUDecomposition solver(n, LUMethod::CHOLESKY);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else if (choice == 5) {
                GaussJacobi solver(n);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else if (choice == 6) {
                GaussSeidel solver(n);
                solver.loadAugmentedMatrix(Afile, bfile);
                result = solver.solve();
            }
            else {
                cout << "Invalid choice\n";
                return 0;
            }

            cout << "\nSolution:\n";
            for (int i = 0; i < (int)result.size(); i++)
                cout << "x" << i + 1 << " = " << result[i] << "\n";

        } catch (exception& e) {
            cout << "Error: " << e.what() << "\n";
        }
    }
    else if (mainChoice == 2) {
        int r, c;
        cout << "\nEnter rows and columns of Matrix A: ";
        cin >> r >> c;

        Matrix A(r, c);
        cout << "Enter elements of Matrix A:\n";
        cin >> A;

        cout << "\nChoose Operation:\n";
        cout << "1. Addition\n";
        cout << "2. Subtraction\n";
        cout << "3. Multiplication\n";
        cout << "4. Transpose\n";
        cout << "5. Determinant\n";
        cout << "6. Inverse\n";
        cout << "7. Check Properties\n";
        cout << "8. Gershgorin Discs\n";
        cout << "Enter choice: ";

        int op;
        cin >> op;

        try {
            if (op == 1 || op == 2 || op == 3) {
                int r2, c2;
                cout << "\nEnter rows and columns of Matrix B: ";
                cin >> r2 >> c2;
                Matrix B(r2, c2);
                cout << "Enter elements of Matrix B:\n";
                cin >> B;

                if      (op == 1) { Matrix res = A + B; cout << "\nA + B:\n" << res; }
                else if (op == 2) { Matrix res = A - B; cout << "\nA - B:\n" << res; }
                else              { Matrix res = A * B; cout << "\nA * B:\n" << res; }
            }
            else if (op == 4) {
                Matrix res = A.transpose();
                cout << "\nTranspose:\n" << res;
            }
            else if (op == 5) {
                if (!A.isSquare()) throw runtime_error("Determinant requires square matrix");
                cout << "\nDeterminant = " << A.determinant() << "\n";
            }
            else if (op == 6) {
                if (!A.isSquare()) throw runtime_error("Inverse requires square matrix");
                Matrix inv = A.inverse();
                cout << "\nInverse:\n" << inv;
            }
            else if (op == 7) {
                cout << "\nMatrix Properties:\n";
                cout << "Square:              " << A.isSquare()             << "\n";
                cout << "Symmetric:           " << A.isSymmetric()          << "\n";
                cout << "Identity:            " << A.isIdentity()           << "\n";
                cout << "Null:                " << A.isNull()               << "\n";
                cout << "Diagonal:            " << A.isDiagonal()           << "\n";
                cout << "Diagonally Dominant: " << A.isDiagonallyDominant() << "\n";
            }
            else if (op == 8) {
                if (!A.isSquare())
                    throw runtime_error("Gershgorin requires a square matrix");
                Gershgorin g(A);
                g.print();
            }
            else {
                cout << "Invalid operation\n";
            }
        } catch (exception& e) {
            cout << "Error: " << e.what() << "\n";
        }
    }
    else if (mainChoice == 3) {

        cout << "\nChoose Interpolation Method:\n";
        cout << "1. Lagrange Interpolation\n";
        cout << "2. Least Squares Approximation\n";
        cout << "Enter choice: ";

        int interpChoice;
        cin >> interpChoice;

        // ── 1. Lagrange ───────────────────────────────────────────────────
        if (interpChoice == 1) {

            cout << "\nLoad data points from:\n";
            cout << "1. File\n";
            cout << "2. Keyboard\n";
            cout << "Enter choice: ";

            int srcChoice;
            cin >> srcChoice;

            int n;
            LagrangeInterpolation* solver = nullptr;

            try {
                if (srcChoice == 1) {
                    string fname;
                    cout << "Enter filename (e.g. points.txt): ";
                    cin >> fname;

                    ifstream peek(fname);
                    if (!peek) {
                        cout << "Error: cannot open " << fname << "\n";
                        return 1;
                    }
                    peek >> n;
                    peek.close();

                    solver = new LagrangeInterpolation(n);
                    ifstream fin(fname);
                    solver->loadDataPoints(fin);
                }
                else {
                    cout << "Enter number of data points: ";
                    cin >> n;

                    solver = new LagrangeInterpolation(n);

                    cout << "Enter data points as  x  y  pairs:\n";
                    for (int i = 0; i < n; ++i) {
                        double xi, yi;
                        cout << "  Point " << i << ":  x = ";
                        cin >> xi;
                        cout << "           y = ";
                        cin >> yi;
                        (*solver)(i, 0) = xi;
                        (*solver)(i, 1) = yi;
                    }
                    cout << "\nData points accepted.\n\n";
                }

                solver->displayDataPoints();

                bool running = true;
                while (running) {
                    cout << "Lagrange Operations:\n";
                    cout << "1. Interpolate at a point\n";
                    cout << "2. Evaluate all basis polynomials L_i(x)\n";
                    cout << "3. Compute derivative P'(x)\n";
                    cout << "4. Print symbolic formula\n";
                    cout << "5. Back to main menu\n";
                    cout << "Enter choice: ";

                    int op;
                    cin >> op;

                    if (op == 1) {
                        double xq;
                        cout << "Enter x to interpolate: ";
                        cin >> xq;
                        double val = solver->interpolate(xq);
                        cout << "\nP(" << xq << ") = " << val << "\n\n";
                    }
                    else if (op == 2) {
                        double xq;
                        cout << "Enter x to evaluate basis polynomials: ";
                        cin >> xq;
                        vector<double> L = solver->basisPolynomials(xq);
                        cout << "\nBasis polynomials at x = " << xq << ":\n";
                        for (int i = 0; i < (int)L.size(); i++)
                            cout << "  L_" << i << "(" << xq << ") = " << L[i] << "\n";
                        cout << "\n";
                    }
                    else if (op == 3) {
                        double xq, h;
                        cout << "Enter x for derivative: ";
                        cin >> xq;
                        cout << "Enter step h (e.g. 1e-5): ";
                        cin >> h;
                        
                    }
                    else if (op == 4) {
                        solver->printFormula();
                    }
                    else if (op == 5) {
                        running = false;
                    }
                    else {
                        cout << "Invalid choice\n";
                    }
                }

                delete solver;
            }
            catch (exception& e) {
                cout << "Error: " << e.what() << "\n";
                delete solver;
            }
        }

        // ── 2. Least Squares ──────────────────────────────────────────────
        else if (interpChoice == 2) {

            cout << "\nLoad data points from:\n";
            cout << "1. File\n";
            cout << "2. Keyboard\n";
            cout << "Enter choice: ";

            int srcChoice;
            cin >> srcChoice;

            int n, degree;
            LeastSquares* solver = nullptr;

            try {
                cout << "Enter polynomial degree (must be < number of points): ";
                cin >> degree;

                if (srcChoice == 1) {
                    string fname;
                    cout << "Enter filename (e.g. points.txt): ";
                    cin >> fname;

                    ifstream peek(fname);
                    if (!peek) { cout << "Error: cannot open " << fname << "\n"; return 1; }
                    peek >> n;
                    peek.close();

                    solver = new LeastSquares(n, degree);
                    ifstream fin(fname);
                    solver->loadDataPoints(fin);   // also calls fit()
                }
                else {
                    cout << "Enter number of data points: ";
                    cin >> n;

                    solver = new LeastSquares(n, degree);

                    cout << "Enter data points as  x  y  pairs:\n";
                    for (int i = 0; i < n; ++i) {
                        double xi, yi;
                        cout << "  Point " << i << ":  x = "; cin >> xi;
                        cout << "           y = "; cin >> yi;
                        (*solver)(i, 0) = xi;
                        (*solver)(i, 1) = yi;
                    }
                    cout << "\nData points accepted.\n\n";
                    solver->fit();   // must call manually when data entered by keyboard
                }

                solver->displayDataPoints();

                bool running = true;
                while (running) {
                    cout << "Least Squares Operations:\n";
                    cout << "1. Evaluate P(x) at a point\n";
                    cout << "2. Print fitted polynomial formula\n";
                    cout << "3. Show normal equations\n";
                    cout << "4. Show residuals\n";
                    cout << "5. Show SSE (Sum of Squared Errors)\n";
                    cout << "6. Back to main menu\n";
                    cout << "Enter choice: ";

                    int op;
                    cin >> op;

                    if (op == 1) {
                        double xq;
                        cout << "Enter x to evaluate: ";
                        cin >> xq;
                        cout << "\nP(" << xq << ") = " << solver->interpolate(xq) << "\n\n";
                    }
                    else if (op == 2) { solver->printFormula(); }
                    else if (op == 3) { solver->printNormalEquations(); }
                    else if (op == 4) {
                        vector<double> res = solver->residuals();
                        cout << "\nResiduals  r_i = y_i - P(x_i):\n";
                        for (int i = 0; i < (int)res.size(); ++i)
                            cout << "  r[" << i << "] = " << res[i] << "\n";
                        cout << "\n";
                    }
                    else if (op == 5) { cout << "\nSSE = " << solver->SSE() << "\n\n"; }
                    else if (op == 6) { running = false; }
                    else { cout << "Invalid choice\n"; }
                }

                delete solver;
            }
            catch (exception& e) {
                cout << "Error: " << e.what() << "\n";
                delete solver;
            }
        }

        else {
            cout << "Invalid interpolation choice\n";
        }
    }

    // ─────────────────────────────────────────────
    //  Unknown
    // ─────────────────────────────────────────────
    else {
        cout << "Invalid main choice\n";
    }

    return 0;
}