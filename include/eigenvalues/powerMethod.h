#ifndef POWER_METHOD_H
#define POWER_METHOD_H

#include "matrixUtils.h"

// Método da Potência Regular
// Encontra o autovalor dominante (maior |lambda|) e o autovetor correspondente.
// Entradas: A (n x n), v0 (vetor inicial), eps (tolerância)
// Saídas: lambda (autovalor dominante), x (autovetor)
void regularPower(const Matrix &A, const Vector &v0, double eps,
                  double &lambda, Vector &x);

// Método da Potência Inverso
// Encontra o autovalor de menor |lambda| e seu autovetor.
// Usa decomposição LU para evitar calcular a inversa.
void inversePower(const Matrix &A, const Vector &v0, double eps,
                  double &lambda, Vector &x);

// Método da Potência com Deslocamento
// Encontra o autovalor mais próximo de mu e seu autovetor.
// Constrói A_hat = A - mu*I e aplica inversePower.
void shiftedPower(const Matrix &A, const Vector &v0, double eps, double mu,
                  double &lambda, Vector &x);

#endif