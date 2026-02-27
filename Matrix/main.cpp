#include <iostream>
#include <fstream>
#include <vector>
#include "../Include/Matrix.hpp"
#include "../Include/GaussianElimination.hpp"

using namespace std;

int main()
{
    // Addition and Subtraction input files
    ifstream matrixInput("input_matrix.txt");
    ofstream matrixOutput("Add_sub_output.txt");

    if (matrixInput)
    {
        int r, c;
        matrixInput >> r >> c;

        Matrix A(r, c);
        Matrix B(r, c);

        // Read Matrix A
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
            {
                double val;
                matrixInput >> val;
                A.setValue(i, j, val);
            }

        // Read Matrix B
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
            {
                double val;
                matrixInput >> val;
                B.setValue(i, j, val);
            }

        Matrix C = A.add(B);
        Matrix D = A.subtract(B);

        matrixOutput << "Addition Result:\n";
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
                matrixOutput << C.getValue(i, j) << " ";
            matrixOutput << endl;
        }

        matrixOutput << "\nSubtraction Result:\n";
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
                matrixOutput << D.getValue(i, j) << " ";
            matrixOutput << endl;
        }

        cout << "Addition and Subtraction stored in Add_sub_output.txt\n";

        matrixInput.close();
        matrixOutput.close();
    }
    else
    {
        cout << "Matrix input file not found.\n";
    }

// Gaussian Elimination input files
    ifstream Afile("49l.txt");   // coefficient matrix (49x49)
    ifstream bfile("49r.txt");   // RHS vector (49x1)
    ofstream gOutput("49_solution.txt");

    if (!Afile || !bfile)
    {
        cout << "Error opening 49l.txt or 49r.txt\n";
        return 1;
    }

    int n, m;
    Afile >> n >> m;   // expected: 49 49

    GaussianElimination G(n, m);

    // Combine A and b into augmented matrix [A | b]
    G.loadAugmentedMatrix(Afile, bfile);

    vector<double> solution = G.solve();

    for (int i = 0; i < n; i++)
        gOutput << "x" << i + 1 << " = " << solution[i] << endl;

    Afile.close();
    bfile.close();
    gOutput.close();

    cout << "Gaussian solution stored in 49_solution.txt\n";

    return 0;
}