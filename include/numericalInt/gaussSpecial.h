#ifndef GAUSSSPECIAL_H
#define GAUSSSPECIAL_H

#include <functional>

using Func1D = std::function<double(double)>;

// QUADRATURAS ESPECIAIS DE GAUSS
// Gauss-Hermite:   I = ∫_{-∞}^{+∞} e^{-s²} f(s) ds ≈ Σ w_k f(s_k)
// Gauss-Laguerre:  I = ∫_{0}^{+∞}  e^{-s}  f(s) ds ≈ Σ w_k f(s_k)
// Gauss-Chebyshev: I = ∫_{-1}^{+1} 1/√(1-s²) f(s) ds ≈ Σ w_k f(s_k)
// n = 2, 3 ou 4  (sem partição)

double gaussHermite(Func1D f, int n);
double gaussLaguerre(Func1D f, int n);
double gaussChebyshev(Func1D f, int n);

#endif