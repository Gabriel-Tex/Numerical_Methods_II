#include "../../include/eigenValues/metodoQR.h"
#include <cmath>

// Matriz de rotacao para decomposicao QR
// tg(theta) = -A[i][j] / A[j][j] 
// Indices i e j sao 0-based.
Matriz matrizJacobiBaseadaNoElemento_ij_DeRvelha(const Matriz &A,
                                                 int i, int j, int n)
{
    const double eps = 1e-9;
    Matriz Gij = matrizIdentidade(n);

    if (std::fabs(A[i][j]) <= eps)
        return Gij;

    double theta;

    if (std::fabs(A[j][j]) <= eps)
    {
        theta = (A[i][j] < 0.0) ? M_PI / 2.0 : -M_PI / 2.0;
    }
    else
    {
        theta = std::atan(-A[i][j] / A[j][j]);
    }

    double c = std::cos(theta);
    double s = std::sin(theta);

    Gij[j][j] = c;
    Gij[i][i] = c;
    Gij[i][j] = s;
    Gij[j][i] = -s;

    return Gij;
}

// Aula 22 — Algoritmo decomposicaoQR (secao 3.1.2)
// A = Q * R, onde Q e ortogonal e R e triangular superior.
// Q^T e montado como produto das rotacoes aplicadas a esquerda.
void decomposicaoQR(const Matriz &A, int n, Matriz &Q, Matriz &R)
{
    Matriz QT = matrizIdentidade(n);
    Matriz R_prev = A;

    for (int j = 0; j <= n - 2; j++)
    {
        for (int i = j + 1; i <= n - 1; i++)
        {
            Matriz Gij = matrizJacobiBaseadaNoElemento_ij_DeRvelha(R_prev, i, j, n);

            // R_nova = Gij * R_prev  (multiplicacao a esquerda)
            R_prev = produtoMatrizes(Gij, R_prev);

            // QT <- Gij * QT
            QT = produtoMatrizes(Gij, QT);
        }
    }

    Q = transposta(QT);
    R = R_prev;
}

// Aula 22 — Algoritmo 3.1.1: metodo QR (loop de iteracoes)
void metodoQR(const Matriz &A,
              int n,
              double eps,
              Matriz &P,
              Vetor &Lamb)
{
    P = matrizIdentidade(n);
    Matriz A_prev = A;
    double val = 100.0;

    while (val > eps)
    {
        Matriz Q, R;
        decomposicaoQR(A_prev, n, Q, R);

        // A_nova = R * Q  (ordem reversa — equivale a Q^T * A_prev * Q)
        Matriz A_nova = produtoMatrizes(R, Q);
        A_prev = A_nova;

        // P <- P * Q
        P = produtoMatrizes(P, Q);

        val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A_nova, n);
    }

    Lamb.resize(n);
    for (int i = 0; i < n; i++)
        Lamb[i] = A_prev[i][i];
}