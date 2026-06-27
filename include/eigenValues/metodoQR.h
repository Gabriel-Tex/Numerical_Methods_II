#ifndef METODO_QR_H
#define METODO_QR_H

#include "matrizUtils.h"

// Método QR
//
// Diagonaliza a matriz simétrica A através de múltiplas
// iterações QR, até que os elementos abaixo da diagonal
// sejam menores que eps.
//
// Cada iteração: decompõe A_prev = Q*R e define A_nova = R*Q.
// As matrizes Q são acumuladas para recuperar autovetores.
//
// Retorna:
//   P - matriz acumulada dos Q's (colunas = autovetores de A)
//   Lamb - vetor com os autovalores de A (diagonal final)
void metodoQR(const Matriz &A,
              int n,
              double eps,
              Matriz &P,
              Vetor &Lamb);

// Decomposição QR de uma matriz A.
// Retorna Q (ortogonal) e R (triangular superior)
// tais que A = Q * R.
//
// Usa rotações de Jacobi para zerar os elementos abaixo
// da diagonal de A (varredura similar à usada no método de Jacobi,
// porém sem a pré-multiplicação pela transposta - é apenas Q^T * A).
void decomposicaoQR(const Matriz &A, int n, Matriz &Q, Matriz &R);

// Constrói a matriz de rotação de Jacobi G_ij para uso na decomposição QR.
// tg(θ) = -A[i][j] / A[j][j]
// Índices i e j são 0-based.
Matriz matrizJacobiBaseadaNoElemento_ij_DeRvelha(const Matriz &A,
                                                 int i, int j, int n);

#endif