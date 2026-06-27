#ifndef SINGULARINTEGRAL_H
#define SINGULARINTEGRAL_H

#include <functional>
#include <string>

using Func1D = std::function<double(double)>;

// INTEGRAÇÃO COM SINGULARIDADE
//
// Estratégia: mudança de variável exponencial (simples ou dupla)
// que mapeia [a,b] → (-∞,+∞), transformando o integrando
// em uma curva em forma de dino adormecido (kkkkkkk).
//
// Solução 2: aproximar a integral infinita por [-c,+c]
// e usar Gauss-Legendre ou Newton-Cotes nesse intervalo,
// aumentando c até convergir.

// Mudança de variável - exponencial simples
// x(s) = (a+b)/2 + (b-a)/2 * tanh(s)
// dx/ds = (b-a)/2 * sech²(s)
double xSimples(double s, double a, double b);
double dxSimples(double s, double a, double b);

// Mudança de variável - exponencial dupla
// x(s) = (a+b)/2 + (b-a)/2 * tanh(π/2 * sinh(s))
// dx/ds = (b-a)/2 * [π/2 * cosh(s) / cosh²(π/2 * sinh(s))]
double xDupla(double s, double a, double b);
double dxDupla(double s, double a, double b);

// Solução 2: integral de f_bar(s) de [-c,+c] usando Gauss-Legendre (n pontos)
// Aumenta c (cMin, cMin+1, ...) até |I_novo - I_ant| < tolerancia
// 'tipo': "simples" ou "dupla"
double singularSolucao2(Func1D f, double a, double b, const std::string &tipo,
                        double cMin, double tolerancia, int nGL, double &cFinal);

#endif