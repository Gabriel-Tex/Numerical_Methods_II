#include "../../include/eigenvalues/qrMethod.h"
#include <cmath>

Matrix givensRotation(const Matrix &A, int i, int j, int n)
{
    const double eps = 1e-9;
    Matrix Gij = identityMatrix(n);

    if (std::fabs(A[i][j]) <= eps)
        return Gij;

    double theta;
    if (std::fabs(A[j][j]) <= eps)
        theta = (A[i][j] < 0.0) ? M_PI / 2.0 : -M_PI / 2.0;
    else
        theta = std::atan(-A[i][j] / A[j][j]);

    double c = std::cos(theta);
    double s = std::sin(theta);

    Gij[j][j] = c;
    Gij[i][i] = c;
    Gij[i][j] = s;
    Gij[j][i] = -s;

    return Gij;
}

void qrDecomposition(const Matrix &A, int n, Matrix &Q, Matrix &R)
{
    Matrix QT = identityMatrix(n);
    Matrix Rprev = A;

    for (int j = 0; j <= n - 2; j++)
    {
        for (int i = j + 1; i <= n - 1; i++)
        {
            Matrix Gij = givensRotation(Rprev, i, j, n);
            Rprev = matMatProduct(Gij, Rprev);
            QT = matMatProduct(Gij, QT);
        }
    }

    Q = transpose(QT);
    R = Rprev;
}

void qrMethod(const Matrix &A, int n, double eps,
              Matrix &P, Vector &eigenvalues)
{
    P = identityMatrix(n);
    Matrix Aprev = A;
    double val = 100.0;

    while (val > eps)
    {
        Matrix Q, R;
        qrDecomposition(Aprev, n, Q, R);
        // A_new = R * Q (equivale a Q^T * A_prev * Q)
        Matrix Anew = matMatProduct(R, Q);
        Aprev = Anew;
        P = matMatProduct(P, Q);
        val = sumSquaresBelowDiagonal(Anew, n);
    }

    eigenvalues.resize(n);
    for (int i = 0; i < n; i++)
        eigenvalues[i] = Aprev[i][i];
}