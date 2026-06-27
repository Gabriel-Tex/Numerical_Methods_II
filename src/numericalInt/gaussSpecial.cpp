#include "../../include/numericalInt/gaussSpecial.h"
#include <cmath>
#include <stdexcept>

// Raízes e pesos - Gauss-Hermite
// n=2: s=±1/√2,  w=√π/2
// n=3: s=0,±√(3/2),  w=2√π/3, √π/6
// n=4: s=±α, ±β
static void hermiteData(int n, double *s, double *w)
{
    switch (n)
    {
    case 2:
        s[0] = -1.0 / sqrt(2.0);
        w[0] = sqrt(M_PI) / 2.0;
        s[1] = +1.0 / sqrt(2.0);
        w[1] = sqrt(M_PI) / 2.0;
        break;
    case 3:
        s[0] = -sqrt(3.0 / 2.0);
        w[0] = sqrt(M_PI) / 6.0;
        s[1] = 0.0;
        w[1] = 2.0 * sqrt(M_PI) / 3.0;
        s[2] = +sqrt(3.0 / 2.0);
        w[2] = sqrt(M_PI) / 6.0;
        break;
    case 4:
        // Raízes de H4(s)=16s⁴-48s²+12=0 → s²=(3±√6)/2 / (raízes do polinômio de Hermite grau 4)
        // s = ±√((3-√6)/2), ±√((3+√6)/2)
        s[0] = -sqrt((3.0 + sqrt(6.0)) / 2.0);
        w[0] = sqrt(M_PI) / (4.0 * (3.0 + sqrt(6.0)));
        s[1] = -sqrt((3.0 - sqrt(6.0)) / 2.0);
        w[1] = sqrt(M_PI) / (4.0 * (3.0 - sqrt(6.0)));
        s[2] = +sqrt((3.0 - sqrt(6.0)) / 2.0);
        w[2] = w[1];
        s[3] = +sqrt((3.0 + sqrt(6.0)) / 2.0);
        w[3] = w[0];
        break;
    default:
        throw std::invalid_argument("Gauss-Hermite: n deve ser 2, 3 ou 4");
    }
}

// Raízes e pesos - Gauss-Laguerre
// n=2: s=2±√2,   w=(2±√2)/4
// n=3: s≈0.4158, 2.2943, 6.2900,  w≈0.7111, 0.2785, 0.0104
// n=4: valores tabelados
static void laguerreData(int n, double *s, double *w)
{
    switch (n)
    {
    case 2:
        s[0] = 2.0 - sqrt(2.0);
        w[0] = (2.0 + sqrt(2.0)) / 4.0;
        s[1] = 2.0 + sqrt(2.0);
        w[1] = (2.0 - sqrt(2.0)) / 4.0;
        break;
    case 3:
        s[0] = 0.4157745568;
        w[0] = 0.7110930099;
        s[1] = 2.2942803603;
        w[1] = 0.2785177336;
        s[2] = 6.2899450829;
        w[2] = 0.0103892565;
        break;
    case 4:
        s[0] = 0.3225476896;
        w[0] = 0.6031541043;
        s[1] = 1.7457611012;
        w[1] = 0.3574186924;
        s[2] = 4.5366202969;
        w[2] = 0.0388879085;
        s[3] = 9.3950709123;
        w[3] = 0.0005392947;
        break;
    default:
        throw std::invalid_argument("Gauss-Laguerre: n deve ser 2, 3 ou 4");
    }
}

// Raízes e pesos - Gauss-Chebyshev
// s_k = cos((2k-1)π / (2n)),  w_k = π/n para todo k
static void chebyshevData(int n, double *s, double *w)
{
    if (n < 2 || n > 4)
        throw std::invalid_argument("Gauss-Chebyshev: n deve ser 2, 3 ou 4");
    for (int k = 1; k <= n; k++)
    {
        s[k - 1] = cos((2.0 * k - 1.0) * M_PI / (2.0 * n));
        w[k - 1] = M_PI / n;
    }
}

// Gauss-Hermite: I = ∫ e^{-s²} f(s) ds ≈ Σ w_k f(s_k)
double gaussHermite(Func1D f, int n)
{
    double s[4], w[4];
    hermiteData(n, s, w);
    double soma = 0.0;
    for (int k = 0; k < n; k++)
        soma += w[k] * f(s[k]);
    return soma;
}

// Gauss-Laguerre: I = ∫ e^{-s} f(s) ds ≈ Σ w_k f(s_k)
double gaussLaguerre(Func1D f, int n)
{
    double s[4], w[4];
    laguerreData(n, s, w);
    double soma = 0.0;
    for (int k = 0; k < n; k++)
        soma += w[k] * f(s[k]);
    return soma;
}

// Gauss-Chebyshev: I = ∫ 1/√(1-s²) f(s) ds ≈ Σ w_k f(s_k)
double gaussChebyshev(Func1D f, int n)
{
    double s[4], w[4];
    chebyshevData(n, s, w);
    double soma = 0.0;
    for (int k = 0; k < n; k++)
        soma += w[k] * f(s[k]);
    return soma;
}