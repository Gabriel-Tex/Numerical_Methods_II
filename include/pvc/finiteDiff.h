#ifndef FINITE_DIFF_H
#define FINITE_DIFF_H

#include <vector>

// PVC1 por Diferenças Finitas com N partições
// PVC1: y'' - y = 0,  y(0)=0, y(1)=1
// Retorna os valores aproximados nos nós interiores x_1..x_{N-1}
std::vector<double> finiteDiffPVC1(int N);

// PVC2 por Diferenças Finitas com N x N partições
// PVC2: d2u/dx2 + d2u/dy2 = 4,  u=0 no contorno
// Retorna os (N-1)^2 valores nos nós interiores, por linha (y crescente)
std::vector<double> finiteDiffPVC2(int N);

#endif