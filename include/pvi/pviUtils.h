#ifndef PVI_UTILS_H
#define PVI_UTILS_H

#include <vector>
#include <functional>

// Estado S(t): vetor de variáveis dependentes
using State = std::vector<double>;

// Tipo da função F(S, t): recebe estado e tempo, retorna dS/dt
using StateFunc = std::function<State(const State &S, double t)>;

// Soma elemento a elemento de dois estados
State addStates(const State &a, const State &b);

// Multiplica estado por escalar
State scaleState(double alpha, const State &a);

// Erro relativo entre dois estados (norma infinito relativa)
double relativeError(const State &newer, const State &older);

#endif