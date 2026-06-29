#ifndef EXPLICIT_EULER_H
#define EXPLICIT_EULER_H

#include "pviUtils.h"
#include <vector>

// Método de Euler Explícito (Forward Euler)
// S_{i+1} = S_i + dt * F(S_i, t_i)
// Retorna todos os estados: states[0]=S0, states[1]=S1, ...
std::vector<State> explicitEuler(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times);

#endif