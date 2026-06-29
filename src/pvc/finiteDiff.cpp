#include "../../include/pvc/finiteDiff.h"
#include "../../include/pvc/pvcUtils.h"

// PVC1: discretização central 
//   (1/dx^2)*y_{i-1} - (2/dx^2+1)*y_i + (1/dx^2)*y_{i+1} = 0
std::vector<double> finiteDiffPVC1(int N)
{
    int M = N - 1;
    double dx = 1.0 / N;
    double dx2 = dx * dx;
    double coefLat = 1.0 / dx2;
    double coefCen = -(2.0 / dx2 + 1.0);

    std::vector<double> a(M, coefLat);
    std::vector<double> b(M, coefCen);
    std::vector<double> c(M, coefLat);
    std::vector<double> d(M, 0.0);
    d[M - 1] = -coefLat * 1.0; // condição y(1) = 1

    return solveTridiagonal(a, b, c, d);
}

// PVC2: máscara de 5 pontos 
//   nó (i,j): index = j*M + i,  M = N-1
std::vector<double> finiteDiffPVC2(int N)
{
    int M = N - 1;
    int total = M * M;
    double dx = 1.0 / N;
    double dx2 = dx * dx;
    double coefLat = 1.0 / dx2;
    double coefCen = -4.0 / dx2;
    double fval = 4.0;

    std::vector<std::vector<double>> A(total, std::vector<double>(total, 0.0));
    std::vector<double> rhs(total, fval);

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < M; i++)
        {
            int k = j * M + i;
            A[k][k] = coefCen;
            if (i > 0)
                A[k][k - 1] = coefLat;
            if (i < M - 1)
                A[k][k + 1] = coefLat;
            if (j > 0)
                A[k][k - M] = coefLat;
            if (j < M - 1)
                A[k][k + M] = coefLat;
        }
    }

    return solveGauss(A, rhs);
}