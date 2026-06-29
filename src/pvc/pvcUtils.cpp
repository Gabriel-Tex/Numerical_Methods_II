#include "../../include/pvc/pvcUtils.h"
#include <cmath>
#include <stdexcept>

std::vector<double> solveTridiagonal(
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> c,
    std::vector<double> d)
{
    int n = (int)b.size();
    std::vector<double> x(n);
    for (int i = 1; i < n; i++)
    {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    return x;
}

std::vector<double> solveGauss(
    std::vector<std::vector<double>> A,
    std::vector<double> b)
{
    int n = (int)b.size();
    for (int col = 0; col < n; col++)
    {
        int pivot = col;
        for (int row = col + 1; row < n; row++)
            if (std::fabs(A[row][col]) > std::fabs(A[pivot][col]))
                pivot = row;
        std::swap(A[col], A[pivot]);
        std::swap(b[col], b[pivot]);
        if (std::fabs(A[col][col]) < 1e-15)
            throw std::runtime_error("Sistema singular");
        double div = A[col][col];
        for (int j = col; j < n; j++)
            A[col][j] /= div;
        b[col] /= div;
        for (int row = 0; row < n; row++)
        {
            if (row == col)
                continue;
            double factor = A[row][col];
            for (int j = col; j < n; j++)
                A[row][j] -= factor * A[col][j];
            b[row] -= factor * b[col];
        }
    }
    return b;
}

double exactSolutionPVC1(double x)
{
    double den = std::exp(-1.0) - std::exp(1.0);
    return (std::exp(-x) - std::exp(x)) / den;
}