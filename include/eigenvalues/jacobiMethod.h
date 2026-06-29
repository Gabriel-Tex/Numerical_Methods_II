#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

#include "matrixUtils.h"

// Diagonaliza a matriz simétrica A por varreduras de Jacobi.
// Retorna P (colunas = autovetores) e eigenvalues (autovalores).
void jacobiMethod(const Matrix &A, int n, double eps,
                  Matrix &P, Vector &eigenvalues);

// Executa uma varredura de Jacobi: percorre todos os pares (i,j)
// com i > j e zera A[i][j] via transformação de similaridade.
// Retorna Anew (matriz após varredura) e Pvar (produto das rotações).
void jacobiSweep(const Matrix &A, int n, Matrix &Anew, Matrix &Pvar);

// Constrói a matriz de rotação de Jacobi G_ij para zerar A[i][j].
// Índices i e j são 0-based.
Matrix jacobiRotation(const Matrix &A, int i, int j, int n);

#endif