#ifndef QR_METHOD_H
#define QR_METHOD_H

#include "matrixUtils.h"

// Diagonaliza a matriz simétrica A por iterações QR.
// Cada iteração: decompõe A_prev = Q*R e define A_new = R*Q.
// Retorna P (colunas = autovetores) e eigenvalues (autovalores).
void qrMethod(const Matrix &A, int n, double eps,
              Matrix &P, Vector &eigenvalues);

// Decomposição QR de A: retorna Q (ortogonal) e R (triangular superior).
// Usa rotações de Givens para zerar os elementos abaixo da diagonal.
void qrDecomposition(const Matrix &A, int n, Matrix &Q, Matrix &R);

// Constrói a matriz de rotação de Givens G_ij para uso na decomposição QR.
// tg(theta) = -A[i][j] / A[j][j]. Índices 0-based.
Matrix givensRotation(const Matrix &A, int i, int j, int n);

#endif