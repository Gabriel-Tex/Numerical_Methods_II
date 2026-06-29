#include "../../include/pvc/finiteElements.h"
#include "../../include/pvc/pvcUtils.h"

// MEF PVC1 
// Matriz local 2x2 
// K_e = (1/h)*[1 -1; -1 1] + (h/6)*[2 1; 1 2]
std::vector<double> finiteElementsPVC1(int N)
{
    int M = N - 1;
    double h = 1.0 / N;
    double K11 = 1.0 / h + 2.0 * h / 6.0;
    double K12 = -1.0 / h + h / 6.0;

    std::vector<double> a(M, K12);
    std::vector<double> b(M, 2.0 * K11);
    std::vector<double> c(M, K12);
    std::vector<double> d(M, 0.0);
    d[M - 1] = -K12 * 1.0; // contribuição de y_N = 1

    return solveTridiagonal(a, b, c, d);
}

// MEF PVC2
// Matriz local 4x4, vetor de carga local fe = -1/16 
// Nós locais: 1=SW, 2=SE, 3=NE, 4=NW
std::vector<double> finiteElementsPVC2(int N)
{
    int M = N - 1;
    int total = M * M;

    double Ke[4][4] = {
        {2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0},
        {-1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0},
        {-1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0},
        {-1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0}};
    double fe = -1.0 / 16.0;

    std::vector<std::vector<double>> A(total, std::vector<double>(total, 0.0));
    std::vector<double> rhs(total, 0.0);

    auto gidx = [&](int i, int j) -> int
    {
        if (i == 0 || i == N || j == 0 || j == N)
            return -1;
        return (j - 1) * M + (i - 1);
    };

    for (int ej = 0; ej < N; ej++)
    {
        for (int ei = 0; ei < N; ei++)
        {
            int gn[4] = {
                gidx(ei, ej),
                gidx(ei + 1, ej),
                gidx(ei + 1, ej + 1),
                gidx(ei, ej + 1)};
            for (int a = 0; a < 4; a++)
            {
                if (gn[a] < 0)
                    continue;
                rhs[gn[a]] += fe;
                for (int b = 0; b < 4; b++)
                {
                    if (gn[b] < 0)
                        continue;
                    A[gn[a]][gn[b]] += Ke[a][b];
                }
            }
        }
    }

    return solveGauss(A, rhs);
}