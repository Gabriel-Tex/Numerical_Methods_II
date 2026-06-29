#ifndef FINITE_ELEMENTS_H
#define FINITE_ELEMENTS_H

#include <vector>

// PVC1 por MEF com N elementos lineares
// PVC1: y'' - y = 0,  y(0)=0, y(1)=1
// Retorna os valores aproximados nos nós interiores 1..N-1
std::vector<double> finiteElementsPVC1(int N);

// PVC2 por MEF com N x N elementos bilineares, f=4
// PVC2: d2u/dx2 + d2u/dy2 = 4,  u=0 no contorno
// Retorna os (N-1)^2 valores nos nós interiores
std::vector<double> finiteElementsPVC2(int N);

#endif