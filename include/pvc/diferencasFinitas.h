#ifndef DIFERENCAS_FINITAS_H
#define DIFERENCAS_FINITAS_H

#include <vector>

// Resolve PVC1 por Diferencas Finitas com N particoes
// PVC1: y'' - y = 0,  y(0) = 0,  y(1) = 1
// Retorna os valores aproximados nos nos interiores x_1, ..., x_{N-1}
std::vector<double> diferencasFinitasPVC1(int N);

// Resolve PVC2 por Diferencas Finitas com N x N particoes
// PVC2: d2u/dx2 + d2u/dy2 = f(x,y),  u = 0 em todo o contorno
// f(x,y) = 4 (constante) em todo o dominio
// Retorna os valores nos (N-1)^2 nos interiores, ordenados linha por linha (y crescente)
std::vector<double> diferencasFinitasPVC2(int N);

#endif