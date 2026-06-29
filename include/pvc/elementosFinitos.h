#ifndef ELEMENTOS_FINITOS_H
#define ELEMENTOS_FINITOS_H

#include <vector>

// Resolve PVC1 por MEF com N elementos (N+1 nos, nos 0 e N sao Dirichlet)
// PVC1: y'' - y = 0,  y(0)=0, y(1)=1
// Retorna os valores aproximados nos nos interiores 1..N-1
std::vector<double> elementosFinitasPVC1(int N);

// Resolve PVC2 por MEF com N x N elementos
// PVC2: d2u/dx2 + d2u/dy2 = f(x,y),  u=0 no contorno, f=4
// Retorna os valores nos (N-1)^2 nos interiores
std::vector<double> elementosFinitasPVC2(int N);

#endif