#ifndef POTENCIA_H
#define POTENCIA_H

#include "matrizUtils.h"

// Método da Potência Regular 
//
// Encontra o autovalor dominante (maior |lambda|) e
// o autovetor correspondente da matriz A.
//
// Entradas: A (n x n), v0 (vetor inicial), eps (tolerância)
// Saída: lambda1 (autovalor dominante), x1 (autovetor)
void potenciaRegular(const Matriz &A,
                     const Vetor &v0,
                     double eps,
                     double &lambda,
                     Vetor &x);

// Método da Potência Inverso 
//
// Encontra o autovalor de menor |lambda| e seu autovetor.
// Usa decomposição LU para evitar calcular a inversa.
//
// Entradas: A, v0, eps
// Saída: lambda_n (menor autovalor), x_n (autovetor)
void potenciaInverso(const Matriz &A,
                     const Vetor &v0,
                     double eps,
                     double &lambda,
                     Vetor &x);

// Método da Potência com Deslocamento 
//
// Encontra o autovalor mais próximo de mu e seu autovetor.
// Constrói A_hat = A - mu*I e aplica potenciaInverso.
//
// Entradas: A, v0, eps, mu (deslocamento)
// Saída: lambda_i (autovalor próximo de mu), x_i (autovetor)
void potenciaComDeslocamento(const Matriz &A,
                             const Vetor &v0,
                             double eps,
                             double mu,
                             double &lambda,
                             Vetor &x);

#endif