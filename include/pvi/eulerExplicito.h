#ifndef EULER_EXPLICITO_H
#define EULER_EXPLICITO_H

#include "pviUtils.h"
#include <vector>

// Metodo de Euler Explicito (Forward Euler)
// Eq. (12): S_{i+1} = S_i + Delta_t * F(S_i, t_i)
//
// Retorna todos os estados calculados: states[0] = S0, states[1] = S1, ...
// Os tempos correspondentes sao armazenados em 'tempos'.
std::vector<Estado> eulerExplicito(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos);

#endif