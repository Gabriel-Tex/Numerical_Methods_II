#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include "pviUtils.h"
#include <vector>

// Resolutor implícito: dado S_i, t_i e dt, retorna S_{i+1}
// Encapsula a solução analítica (ou numérica) do sistema implícito
using ImplicitSolver = std::function<State(const State &Si, double ti, double dt)>;

// Método de Euler Implícito (Backward Euler)
// S_{i+1} = S_i + dt * F(S_{i+1}, t_{i+1})  (resolvido pelo solver)
std::vector<State> implicitEuler(
    const ImplicitSolver &solver,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times);

#endif