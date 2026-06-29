#include "../../include/eigenvalues/jacobiMethod.h"
#include <cmath>

Matrix jacobiRotation(const Matrix &A, int i, int j, int n)
{
    const double eps = 1e-9;
    Matrix Gij = identityMatrix(n);

    if (std::fabs(A[i][j]) <= eps)
        return Gij;

    double theta;
    if (std::fabs(A[i][i] - A[j][j]) <= eps)
        theta = M_PI / 4.0;
    else
        theta = 0.5 * std::atan(-2.0 * A[i][j] / (A[i][i] - A[j][j]));

    double c = std::cos(theta);
    double s = std::sin(theta);

    Gij[j][j] = c;
    Gij[i][i] = c;
    Gij[i][j] = s;
    Gij[j][i] = -s;

    return Gij;
}

void jacobiSweep(const Matrix &A, int n, Matrix &Anew, Matrix &Pvar)
{
    Pvar = identityMatrix(n);
    Matrix Aprev = A;

    for (int j = 0; j <= n - 2; j++)
    {
        for (int i = j + 1; i <= n - 1; i++)
        {
            Matrix Gij = jacobiRotation(Aprev, i, j, n);
            Matrix GijT = transpose(Gij);
            Anew = matMatProduct(matMatProduct(GijT, Aprev), Gij);
            Aprev = Anew;
            Pvar = matMatProduct(Pvar, Gij);
        }
    }

    Anew = Aprev;
}

void jacobiMethod(const Matrix &A, int n, double eps,
                  Matrix &P, Vector &eigenvalues)
{
    P = identityMatrix(n);
    Matrix Aprev = A;
    double val = 100.0;

    while (val > eps)
    {
        Matrix Anew, Pvar;
        jacobiSweep(Aprev, n, Anew, Pvar);
        Aprev = Anew;
        P = matMatProduct(P, Pvar);
        val = sumSquaresBelowDiagonal(Anew, n);
    }

    eigenvalues.resize(n);
    for (int i = 0; i < n; i++)
        eigenvalues[i] = Aprev[i][i];
}