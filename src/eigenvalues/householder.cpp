#include "../../include/eigenvalues/householder.h"
#include <cmath>

Matrix householderMatrix(const Matrix &A, int i)
{
    int n = (int)A.size();

    // Copia os elementos abaixo da diagonal da coluna i (1-based) para w
    Vector w(n, 0.0);
    for (int k = i; k < n; k++)
        w[k] = A[k][i - 1];

    double Lw = vecNorm(w);

    // w' tem Lw na posição i (0-based)
    Vector wPrime(n, 0.0);
    wPrime[i] = Lw;

    // N = w - w'
    Vector N(n, 0.0);
    for (int k = 0; k < n; k++)
        N[k] = w[k] - wPrime[k];

    double LN = vecNorm(N);
    if (LN < 1e-14)
        return identityMatrix(n);

    // n_vec = N / ||N||
    Vector nVec(n);
    for (int k = 0; k < n; k++)
        nVec[k] = N[k] / LN;

    // H = I - 2 * nVec * nVec^T
    Matrix Hi = identityMatrix(n);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            Hi[r][c] -= 2.0 * nVec[r] * nVec[c];

    return Hi;
}

void householderMethod(const Matrix &A, Matrix &Abar, Matrix &H)
{
    int n = (int)A.size();
    H = identityMatrix(n);
    Matrix Aprev = A;

    for (int i = 1; i <= n - 2; i++)
    {
        Matrix Hi = householderMatrix(Aprev, i);
        // Transformação de similaridade: A_i = Hi * A_{i-1} * Hi
        Matrix Ai = matMatProduct(matMatProduct(Hi, Aprev), Hi);
        Aprev = Ai;
        H = matMatProduct(H, Hi);
    }

    Abar = Aprev;
}