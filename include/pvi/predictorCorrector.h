#ifndef PREDICTOR_CORRECTOR_H
#define PREDICTOR_CORRECTOR_H

#include "pviUtils.h"
#include <vector>

// Adams-Bashforth de 2ª ordem (predição-correção)
// Inicialização: S1 por RK2
// Predição:  S_bar = S_i + (dt/2)*(-F_{i-1} + 3*F_i)
// Correção:  S_{i+1} = S_i + (dt/2)*(F_i + F(S_bar, t_{i+1}))
std::vector<State> adamsBashforth2(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    double eps,
    int maxIter,
    std::vector<double> &times);

// Adams-Bashforth-Moulton de 3ª ordem
// Inicialização: S1, S2 por RK3
// Predição:  S_bar = S_i + (dt/12)*(5*F_{i-2} - 16*F_{i-1} + 23*F_i)
// Correção:  S_{i+1} = S_i + (dt/12)*(-F_{i-1} + 8*F_i + 5*F_bar)
std::vector<State> adamsBashforthMoulton3(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    double eps,
    int maxIter,
    std::vector<double> &times);

// Adams-Bashforth-Moulton de 4ª ordem
// Inicialização: S1, S2, S3 por RK4
// Predição:  S_bar = S_i + (dt/24)*(-9*F_{i-3} + 37*F_{i-2} - 59*F_{i-1} + 55*F_i)
// Correção:  S_{i+1} = S_i + (dt/24)*(F_{i-2} - 5*F_{i-1} + 19*F_i + 9*F_bar)
std::vector<State> adamsBashforthMoulton4(
    const StateFunc &F,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    double eps,
    int maxIter,
    std::vector<double> &times);

#endif