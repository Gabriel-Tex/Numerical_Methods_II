#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "matrixUtils.h"

// Transforma a matriz simétrica A em uma matriz tridiagonal A_bar
// via transformações de similaridade de Householder.
// Retorna A_bar e H = H1 * H2 * ... * H_{n-2} (produto acumulado).
// Os autovetores de A se obtêm por: x = H * v_bar
void householderMethod(const Matrix &A, Matrix &Abar, Matrix &H);

// Constrói a matriz de Householder H_i baseada na coluna i de A_{i-1}.
Matrix householderMatrix(const Matrix &A, int i);

#endif