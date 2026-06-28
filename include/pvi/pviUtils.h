#ifndef PVI_UTILS_H
#define PVI_UTILS_H

#include <vector>
#include <functional>

// Estado S(t)
using Estado = std::vector<double>;

// Tipo da função F(S, t): recebe estado e tempo, retorna dS/dt
using FuncaoF = std::function<Estado(const Estado &S, double t)>;

// Soma de dois estados
Estado somaEstado(const Estado &a, const Estado &b);

// Multiplica estado por escalar
Estado escalaEstado(double alpha, const Estado &a);

// Diferença relativa entre dois estados (norma infinito relativa)
double erroRelativo(const Estado &novo, const Estado &antigo);

#endif