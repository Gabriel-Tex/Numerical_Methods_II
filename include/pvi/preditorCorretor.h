#ifndef PREDITOR_CORRETOR_H
#define PREDITOR_CORRETOR_H

#include "pviUtils.h"
#include <vector>

// Metodo Preditor-Corretor de Adams-Bashforth de 2a ordem
//
// Inicializacao: S1 por RK2
// Predicao  [Eq. 19]: S_bar_{i+1} = S_i + (dt/2)*(-F_{i-1} + 3*F_i)
// Correcao  [Eq. 20]: S_{i+1} = S_i + (dt/2)*(F_i + F_{i+1,bar})
//
// Se maxIter > 1, a correcao e repetida ate convergir.
std::vector<Estado> adamsBashforth2(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos);

// Metodo Preditor-Corretor de Adams-Bashforth-Moulton de 3a ordem
//
// Inicializacao: S1, S2 por RK3
// Predicao  [Eq. 43]: S_bar_{i+1} = S_i + (dt/12)*(5*F_{i-2} - 16*F_{i-1} + 23*F_i)
// Correcao  [Eq. 44]: S_{i+1} = S_i + (dt/12)*(-F_{i-1} + 8*F_i + 5*F_{i+1,bar})
std::vector<Estado> adamsBashforthMoulton3(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos);

// Metodo Preditor-Corretor de Adams-Bashforth-Moulton de 4a ordem
//
// Inicializacao: S1, S2, S3 por RK4
// Predicao (derivada a partir do polinomio de Newton de grau 3):
//   S_bar_{i+1} = S_i + (dt/24)*(-9*F_{i-3} + 37*F_{i-2} - 59*F_{i-1} + 55*F_i)
// Correcao:
//   S_{i+1} = S_i + (dt/24)*(F_{i-2} - 5*F_{i-1} + 19*F_i + 9*F_{i+1,bar})
std::vector<Estado> adamsBashforthMoulton4(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos);

#endif