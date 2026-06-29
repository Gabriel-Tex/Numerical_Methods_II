#include "../../include/pvc/pvcUtils.h"
#include <cmath>
#include <stdexcept>
#include <cstdlib>

// Algoritmo de Thomas para sistema tridiagonal
std::vector<double> resolveTridiagonal(
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> c,
    std::vector<double> d)
{
    int n = (int)b.size();
    std::vector<double> x(n);

    // Eliminacao direta
    for (int i = 1; i < n; i++)
    {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    // Substituicao retroativa
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];

    return x;
}

// Eliminacao de Gauss com pivotamento parcial
std::vector<double> resolveGauss(
    std::vector<std::vector<double>> A,
    std::vector<double> b)
{
    int n = (int)b.size();

    for (int col = 0; col < n; col++)
    {
        // Pivotamento
        int pivo = col;
        for (int row = col + 1; row < n; row++)
            if (std::fabs(A[row][col]) > std::fabs(A[pivo][col]))
                pivo = row;

        std::swap(A[col], A[pivo]);
        std::swap(b[col], b[pivo]);

        if (std::fabs(A[col][col]) < 1e-15)
            throw std::runtime_error("Sistema singular");

        double div = A[col][col];
        for (int j = col; j < n; j++) A[col][j] /= div;
        b[col] /= div;

        for (int row = 0; row < n; row++)
        {
            if (row == col) continue;
            double fator = A[row][col];
            for (int j = col; j < n; j++)
                A[row][j] -= fator * A[col][j];
            b[row] -= fator * b[col];
        }
    }

    return b;
}

// Solucao exata: y(x) = 1/(e^{-1} - e) * (e^{-x} - e^x)
double solucaoExataPVC1(double x)
{
    double den = std::exp(-1.0) - std::exp(1.0);
    return (1.0 / den) * (std::exp(-x) - std::exp(x));
}