#ifndef GAUSSLEGENDRE_H
#define GAUSSLEGENDRE_H

#include <functional>

using Func1D = std::function<double(double)>;
using Func2D = std::function<double(double, double)>;

// GAUSS-LEGENDRE
// Quadratura sobre intervalo [xi, xf] com n=2, 3 ou 4 pontos.
// Mudança de variável: x(α) = (xi+xf)/2 + (xf-xi)/2 * α
// Sem partição: função pura de n pontos.
// Com partição: itera N=1,2,4,... até tolerância.

// Uma partição, n pontos (n = 2, 3 ou 4)
double gaussLegendre(Func1D f, double xi, double xf, int n);

// Com particionamento adaptativo
double gaussLegendreParticionado(Func1D f, double a, double b,
                                 int n, double tolerancia, int &iteracoes);

// Versão 2D para integrais duplas sobre domínio retangular [xa,xb] x [ya,yb]
// Usa n pontos em cada direção
double gaussLegendre2D(Func2D f, double xa, double xb, double ya, double yb, int n);

#endif