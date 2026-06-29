#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "pviUtils.h"
#include <vector>

// Runge-Kutta de 2ª ordem
// S_bar = S_i + dt * F(S_i, t_i)
// S_{i+1} = S_i + dt * (1/2*F1 + 1/2*F2)
std::vector<State> rungeKutta2(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times);

// Runge-Kutta de 3ª ordem
std::vector<State> rungeKutta3(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times);

// Runge-Kutta de 4ª ordem
std::vector<State> rungeKutta4(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times);

#endif