#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "matrizUtils.h"

// Método de Householder
//
// Transforma a matriz simétrica A em uma matriz tridiagonal
// A_barra usando transformações de similaridade sucessivas.
//
// Retorna:
//   A_barra - matriz tridiagonal resultante
//   H - produto acumulado das matrizes de Householder
//      (H = I * H1 * H2 * ... * H_{n-2})
//       Usado para recuperar os autovetores de A a
//       partir dos autovetores de A_barra: x = H * v
void metodoDeHouseholder(const Matriz &A,
                         Matriz &A_barra,
                         Matriz &H);

// Constrói a matriz de Householder H_i baseada na coluna i
// da matriz do passo anterior (A_{i-1}).
// A estrutura é: H_i = [ I  0 ]
//                      [ 0  H_II ]
// para preservar as primeiras (i-1) linhas e colunas.
Matriz matrizHouseholderBaseadaNaCol(const Matriz &A, int i);

#endif