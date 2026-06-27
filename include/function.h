#ifndef FUNCTION_H
#define FUNCTION_H

#include <cmath>

// F(x) = x^3
double F(double x);

// Aulas 7, 8 e 11 — função do exemplo de Newton-Cotes e Gauss-Legendre
// f(x) = (sin(2x) + 4x^2 + 3x)^2
double F_NC(double x);

// Aula 14 — funções com singularidade
// f(x) = x^(-1/3),  integral de 0 a 1 = 3/2
double F_sing1(double x);

// f(x) = x^(-1/2) * e^(-x),  integral de 0 a +inf = sqrt(pi)
double F_sing2(double x);

// Aulas 15 e 16 — paraboloide hiperbólico z = 0.2*(x^2 - y^2)
double F_surf(double x, double y);

// Derivadas parciais de F_surf (necessárias para área de superfície 3D)
double dF_surf_dx(double x, double y);
double dF_surf_dy(double x, double y);

#endif