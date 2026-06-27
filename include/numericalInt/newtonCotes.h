#ifndef NEWTONCOTES_H
#define NEWTONCOTES_H

#include <functional>

// NEWTON-COTES
// Cada função de quadratura recebe o integrando como std::function,
// os limites xi e xf, e retorna a estimativa para UMA partição.
// As funções "partitioned" subdividem o intervalo em N partes e
// somam os resultados, iterando até atingir a tolerância.

using Func1D = std::function<double(double)>;

// ---- FECHADAS ----
double trapezioFechado(Func1D f, double xi, double xf);
double simpson13(Func1D f, double xi, double xf);
double simpson38(Func1D f, double xi, double xf);
double fechadaGrau4(Func1D f, double xi, double xf);

// ---- ABERTAS ----
double trapezioAberto(Func1D f, double xi, double xf);
double milne(Func1D f, double xi, double xf);
double abertaGrau3(Func1D f, double xi, double xf);
double abertaGrau4(Func1D f, double xi, double xf);

// ---- COM PARTICIONAMENTO ----
// Itera N=1,2,4,8,... até |I_novo - I_anterior| < tolerancia
// Retorna o resultado e preenche 'iteracoes' com o N final usado
double newtonCotesParticionado(Func1D f, double a, double b,
                               int grau, bool fechada,
                               double tolerancia, int &iteracoes);

#endif