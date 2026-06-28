#ifndef EULER_IMPLICITO_H
#define EULER_IMPLICITO_H

#include "pviUtils.h"
#include <vector>

// Metodo de Euler Implicito (Backward Euler)
// Eq. (17): S_{i+1} = S_i + Delta_t * F(S_{i+1}, t_{i+1})
//
// Para o PVI-2 (sistema linear), a formula pode ser resolvida explicitamente.
// A funcao resolutora recebe S_i, t_i, dt e retorna S_{i+1}.
using ResolutorImplicito = std::function<Estado(const Estado &Si, double ti, double dt)>;

// Resolve usando um resolutor implicito especifico para o PVI em questao.
std::vector<Estado> eulerImplicito(
    const ResolutorImplicito &resolver,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos);

#endif