#include "../../include/eigenvalues/matrixUtils.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

Matrix identityMatrix(int n)
{
    Matrix I(n, Vector(n, 0.0));
    for (int i = 0; i < n; i++)
        I[i][i] = 1.0;
    return I;
}

Matrix transpose(const Matrix &A)
{
    int n = (int)A.size();
    int m = (int)A[0].size();
    Matrix T(m, Vector(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            T[j][i] = A[i][j];
    return T;
}

Vector matVecProduct(const Matrix &A, const Vector &v)
{
    int n = (int)A.size();
    Vector result(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i] += A[i][j] * v[j];
    return result;
}

Matrix matMatProduct(const Matrix &A, const Matrix &B)
{
    int n = (int)A.size();
    int m = (int)B[0].size();
    int k = (int)B.size();
    Matrix C(n, Vector(m, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int p = 0; p < k; p++)
                C[i][j] += A[i][p] * B[p][j];
    return C;
}

double vecNorm(const Vector &v)
{
    double s = 0.0;
    for (double x : v)
        s += x * x;
    return std::sqrt(s);
}

double dotProduct(const Vector &a, const Vector &b)
{
    double s = 0.0;
    for (int i = 0; i < (int)a.size(); i++)
        s += a[i] * b[i];
    return s;
}

Vector normalizeVec(const Vector &v)
{
    double norm = vecNorm(v);
    Vector u(v.size());
    for (int i = 0; i < (int)v.size(); i++)
        u[i] = v[i] / norm;
    return u;
}

double sumSquaresBelowDiagonal(const Matrix &A, int n)
{
    double sum = 0.0;
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            sum += A[i][j] * A[i][j];
    return sum;
}

void printMatrix(const Matrix &A, const std::string &name)
{
    int n = (int)A.size();
    if (!name.empty())
        std::cout << name << ":\n";
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

void printVector(const Vector &v, const std::string &name)
{
    if (!name.empty())
        std::cout << name << ": ";
    std::cout << "[";
    for (int i = 0; i < (int)v.size(); i++)
    {
        std::cout << std::fixed << std::setprecision(6) << v[i];
        if (i + 1 < (int)v.size())
            std::cout << "  ";
    }
    std::cout << "]\n";
}

void luDecomposition(const Matrix &A, Matrix &L, Matrix &U)
{
    int n = (int)A.size();
    L = identityMatrix(n);
    U = A;
    for (int k = 0; k < n; k++)
    {
        if (std::fabs(U[k][k]) < 1e-14)
            throw std::runtime_error("luDecomposition: pivo zero — matriz singular");
        for (int i = k + 1; i < n; i++)
        {
            double m = U[i][k] / U[k][k];
            L[i][k] = m;
            for (int j = k; j < n; j++)
                U[i][j] -= m * U[k][j];
        }
    }
}

Vector solveLower(const Matrix &L, const Vector &b)
{
    int n = (int)b.size();
    Vector y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
    }
    return y;
}

Vector solveUpper(const Matrix &U, const Vector &y)
{
    int n = (int)y.size();
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }
    return x;
}

Vector solveLU(const Matrix &L, const Matrix &U, const Vector &b)
{
    return solveUpper(U, solveLower(L, b));
}

Matrix invertMatrix(const Matrix &A)
{
    int n = (int)A.size();
    Matrix L, U;
    luDecomposition(A, L, U);
    Matrix inv(n, Vector(n, 0.0));
    for (int j = 0; j < n; j++)
    {
        Vector ej(n, 0.0);
        ej[j] = 1.0;
        Vector col = solveLU(L, U, ej);
        for (int i = 0; i < n; i++)
            inv[i][j] = col[i];
    }
    return inv;
}