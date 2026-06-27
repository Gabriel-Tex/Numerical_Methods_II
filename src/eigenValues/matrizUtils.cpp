#include "../../include/eigenValues/matrizUtils.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

Matriz matrizIdentidade(int n)
{
    Matriz I(n, Vetor(n, 0.0));
    for (int i = 0; i < n; i++)
        I[i][i] = 1.0;
    return I;
}

Matriz transposta(const Matriz &A)
{
    int n = (int)A.size();
    int m = (int)A[0].size();
    Matriz T(m, Vetor(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            T[j][i] = A[i][j];
    return T;
}

Vetor produtoMatrizVetor(const Matriz &A, const Vetor &v)
{
    int n = (int)A.size();
    Vetor resultado(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            resultado[i] += A[i][j] * v[j];
    return resultado;
}

Matriz produtoMatrizes(const Matriz &A, const Matriz &B)
{
    int n = (int)A.size();
    int m = (int)B[0].size();
    int k = (int)B.size();
    Matriz C(n, Vetor(m, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int p = 0; p < k; p++)
                C[i][j] += A[i][p] * B[p][j];
    return C;
}

double normaVetor(const Vetor &v)
{
    double s = 0.0;
    for (double x : v)
        s += x * x;
    return std::sqrt(s);
}

double produtoInterno(const Vetor &a, const Vetor &b)
{
    double s = 0.0;
    for (int i = 0; i < (int)a.size(); i++)
        s += a[i] * b[i];
    return s;
}

Vetor normalizaVetor(const Vetor &v)
{
    double norm = normaVetor(v);
    Vetor u(v.size());
    for (int i = 0; i < (int)v.size(); i++)
        u[i] = v[i] / norm;
    return u;
}

double somaDosQuadradosDosTermosAbaixoDaDiagonal(const Matriz &A, int n)
{
    double soma = 0.0;
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            soma += A[i][j] * A[i][j];
    return soma;
}

void imprimirMatriz(const Matriz &A, const std::string &nome)
{
    int n = (int)A.size();
    if (!nome.empty())
        std::cout << nome << ":\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << "  [";
        for (int j = 0; j < (int)A[i].size(); j++)
        {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << A[i][j];
            if (j + 1 < (int)A[i].size())
                std::cout << "  ";
        }
        std::cout << " ]\n";
    }
}

void imprimirVetor(const Vetor &v, const std::string &nome)
{
    if (!nome.empty())
        std::cout << nome << ": ";
    std::cout << "[";
    for (int i = 0; i < (int)v.size(); i++)
    {
        std::cout << std::fixed << std::setprecision(6) << v[i];
        if (i + 1 < (int)v.size())
            std::cout << "  ";
    }
    std::cout << "]\n";
}

void decompLU(const Matriz &A, Matriz &L, Matriz &U)
{
    int n = (int)A.size();
    L = matrizIdentidade(n);
    U = A;

    for (int k = 0; k < n; k++)
    {
        if (std::fabs(U[k][k]) < 1e-14)
            throw std::runtime_error("decompLU: pivô zero — matriz singular ou mal condicionada");
        for (int i = k + 1; i < n; i++)
        {
            double m = U[i][k] / U[k][k];
            L[i][k] = m;
            for (int j = k; j < n; j++)
                U[i][j] -= m * U[k][j];
        }
    }
}

Vetor resolverL(const Matriz &L, const Vetor &b)
{
    int n = (int)b.size();
    Vetor y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
    }
    return y;
}

Vetor resolverU(const Matriz &U, const Vetor &y)
{
    int n = (int)y.size();
    Vetor x(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }
    return x;
}

Vetor solverLU(const Matriz &L, const Matriz &U, const Vetor &b)
{
    Vetor y = resolverL(L, b);
    return resolverU(U, y);
}

Matriz calcularInversa(const Matriz &A)
{
    int n = (int)A.size();
    Matriz L, U;
    decompLU(A, L, U);

    Matriz inv(n, Vetor(n, 0.0));
    for (int j = 0; j < n; j++)
    {
        Vetor ej(n, 0.0);
        ej[j] = 1.0;
        Vetor col = solverLU(L, U, ej);
        for (int i = 0; i < n; i++)
            inv[i][j] = col[i];
    }
    return inv;
}