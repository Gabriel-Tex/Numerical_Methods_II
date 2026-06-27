#include "../../include/numericalInt/gaussLegendre.h"
#include <cmath>
#include <stdexcept>
#include <array>

// Raízes (α_k) e pesos (w_k) dos polinômios de Legendre
// n=2: P2(α)=0  →  α = ±1/√3,  w = 1
// n=3: P3(α)=0  →  α = 0, ±√(3/5),  w = 8/9, 5/9
// n=4: P4(α)=0  →  α = ±√(3/7 ∓ 2/7*√(6/5)),  pesos derivados

struct GLData
{
    int n;
    double alpha[4];
    double w[4];
};

static const GLData GL2 = {
    2,
    {-1.0 / sqrt(3.0), +1.0 / sqrt(3.0), 0.0, 0.0},
    {1.0, 1.0, 0.0, 0.0}};

static const GLData GL3 = {
    3,
    {-sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0), 0.0},
    {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 0.0}};

// n=4: raízes de P4(α) = (35α⁴ - 30α² + 3)/8 = 0
// α = ±sqrt((3 ± 2*sqrt(6/5)) / 7)
static const GLData GL4 = {
    4,
    {-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),
     -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
     +sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
     +sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))},
    {(18.0 - sqrt(30.0)) / 36.0,
     (18.0 + sqrt(30.0)) / 36.0,
     (18.0 + sqrt(30.0)) / 36.0,
     (18.0 - sqrt(30.0)) / 36.0}};

static const GLData &getGL(int n)
{
    switch (n)
    {
    case 2:
        return GL2;
    case 3:
        return GL3;
    case 4:
        return GL4;
    default:
        throw std::invalid_argument("Gauss-Legendre: n deve ser 2, 3 ou 4");
    }
}

// Quadratura de Gauss-Legendre em [xi, xf] com n pontos
// I ≈ (xf-xi)/2 * Σ w_k * f(x(α_k))
// onde x(α) = (xi+xf)/2 + (xf-xi)/2 * α
double gaussLegendre(Func1D f, double xi, double xf, int n)
{
    const GLData &gl = getGL(n);
    double mid = (xi + xf) / 2.0;
    double half = (xf - xi) / 2.0;

    double soma = 0.0;
    for (int k = 0; k < n; k++)
    {
        double x = mid + half * gl.alpha[k];
        soma += gl.w[k] * f(x);
    }
    return half * soma;
}

// Particionamento adaptativo
double gaussLegendreParticionado(Func1D f, double a, double b,
                                 int n, double tolerancia, int &iteracoes)
{
    int N = 1;
    double I_ant = gaussLegendre(f, a, b, n);

    for (int iter = 1; iter <= 100000; iter++)
    {
        N *= 2;
        double dx = (b - a) / N;
        double I_novo = 0.0;
        for (int i = 0; i < N; i++)
        {
            double xi = a + i * dx;
            double xf = xi + dx;
            I_novo += gaussLegendre(f, xi, xf, n);
        }
        iteracoes = N;
        if (fabs(I_novo - I_ant) < tolerancia)
        {
            return I_novo;
        }
        I_ant = I_novo;
    }
    return I_ant;
}

// Gauss-Legendre 2D — domínio retangular [xa,xb] x [ya,yb]
// I ≈ |J| * Σ_i Σ_j w_i * w_j * f(x(α_i), y(β_j))
// |J| = (xb-xa)/2 * (yb-ya)/2
double gaussLegendre2D(Func2D f,
                       double xa, double xb,
                       double ya, double yb,
                       int n)
{
    const GLData &gl = getGL(n);
    double midX = (xa + xb) / 2.0;
    double halfX = (xb - xa) / 2.0;
    double midY = (ya + yb) / 2.0;
    double halfY = (yb - ya) / 2.0;
    double jac = halfX * halfY; // |J| = det([[halfX,0],[0,halfY]])

    double soma = 0.0;
    for (int i = 0; i < n; i++)
    {
        double x = midX + halfX * gl.alpha[i];
        for (int j = 0; j < n; j++)
        {
            double y = midY + halfY * gl.alpha[j];
            soma += gl.w[i] * gl.w[j] * f(x, y);
        }
    }
    return jac * soma;
}