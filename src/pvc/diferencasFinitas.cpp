#include "../../include/pvc/diferencasFinitas.h"
#include "../../include/pvc/pvcUtils.h"
#include <cmath>

// ------------------------------------------------------------------
// PVC1: y'' - y = 0,  y(0)=0, y(1)=1
//
// Discretizacao central (Aula 27, eq. 8):
//   (1/dx^2)*y_{i-1} - (2/dx^2 + 1)*y_i + (1/dx^2)*y_{i+1} = 0
//
// com dx = 1/N, nos interiores i = 1..N-1
// ------------------------------------------------------------------
std::vector<double> diferencasFinitasPVC1(int N)
{
    int M = N - 1; // numero de nos interiores
    double dx = 1.0 / N;
    double dx2 = dx * dx;

    double coefLat = 1.0 / dx2;
    double coefCen = -(2.0 / dx2 + 1.0);

    // Sistema tridiagonal
    std::vector<double> a(M, coefLat);  // subdiagonal
    std::vector<double> b(M, coefCen); // diagonal
    std::vector<double> c(M, coefLat); // superdiagonal
    std::vector<double> d(M, 0.0);    // lado direito

    // Condicao de contorno y(0) = 0 nao contribui (coef = 0)
    // Condicao de contorno y(1) = 1 contribui para o ultimo no
    d[M - 1] = -coefLat * 1.0;

    return resolveTridiagonal(a, b, c, d);
}

// ------------------------------------------------------------------
// PVC2: d2u/dx2 + d2u/dy2 = f(x,y),  u=0 no contorno
// f(x,y) = 4
//
// Grade: N x N particoes => (N-1)^2 nos interiores
// Numeracao: no (i,j) com i=coluna (x), j=linha (y), index = j*(N-1)+i
//
// Mascara (Aula 27, eq. 17):
//   (1/dx^2)*(u_esq + u_dir) + (1/dy^2)*(u_inf + u_sup)
//   - 2*(1/dx^2 + 1/dy^2)*u_ij = f_ij
//
// Como dx = dy = 1/N:  coef = N^2
//   N^2*(u_esq+u_dir+u_inf+u_sup) - 4*N^2*u_ij = f_ij
// ------------------------------------------------------------------
std::vector<double> diferencasFinitasPVC2(int N)
{
    int M = N - 1;   // nos interiores por direcao
    int total = M * M;
    double dx = 1.0 / N;
    double dx2 = dx * dx;

    double coefLat = 1.0 / dx2;
    double coefCen = -2.0 * (1.0 / dx2 + 1.0 / dx2); // = -4/dx^2
    double fval = 4.0;

    // Monta sistema geral M^2 x M^2
    std::vector<std::vector<double>> A(total, std::vector<double>(total, 0.0));
    std::vector<double> rhs(total, fval);

    auto idx = [&](int i, int j) { return j * M + i; };

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < M; i++)
        {
            int k = idx(i, j);
            A[k][k] = coefCen;

            // Vizinho esquerdo (i-1,j)
            if (i > 0) A[k][idx(i-1, j)] = coefLat;
            // Vizinho direito (i+1,j)
            if (i < M-1) A[k][idx(i+1, j)] = coefLat;
            // Vizinho inferior (i,j-1)
            if (j > 0) A[k][idx(i, j-1)] = coefLat;
            // Vizinho superior (i,j+1)
            if (j < M-1) A[k][idx(i, j+1)] = coefLat;
            // Contornos (u=0) nao contribuem para rhs
        }
    }

    return resolveGauss(A, rhs);
}