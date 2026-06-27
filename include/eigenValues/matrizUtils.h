#ifndef MATRIZ_UTILS_H
#define MATRIZ_UTILS_H

#include <vector>
#include <string>

// Tipo matriz: vetor de vetores de double (indexado de 0)
using Matriz = std::vector<std::vector<double>>;
using Vetor  = std::vector<double>;

// ---- Operações básicas ----

// Retorna a matriz identidade n x n
Matriz matrizIdentidade(int n);

// Retorna a transposta de A
Matriz transposta(const Matriz &A);

// Produto matriz-vetor: A * v
Vetor produtoMatrizVetor(const Matriz &A, const Vetor &v);

// Produto de duas matrizes: A * B
Matriz produtoMatrizes(const Matriz &A, const Matriz &B);

// Norma euclidiana de um vetor
double normaVetor(const Vetor &v);

// Produto interno de dois vetores
double produtoInterno(const Vetor &a, const Vetor &b);

// Normaliza um vetor (divide pela norma)
Vetor normalizaVetor(const Vetor &v);

// Soma dos quadrados dos elementos abaixo da diagonal principal
double somaDosQuadradosDosTermosAbaixoDaDiagonal(const Matriz &A, int n);

// Imprime uma matriz formatada
void imprimirMatriz(const Matriz &A, const std::string &nome = "");

// Imprime um vetor formatado
void imprimirVetor(const Vetor &v, const std::string &nome = "");

// Decomposição LU (sem pivotamento parcial)
// Retorna L e U tais que A = L*U
void decompLU(const Matriz &A, Matriz &L, Matriz &U);

// Resolve o sistema triangular inferior L*y = b
Vetor resolverL(const Matriz &L, const Vetor &b);

// Resolve o sistema triangular superior U*x = y
Vetor resolverU(const Matriz &U, const Vetor &y);

// Resolve A*x = b usando decomposição LU (L e U pré-calculados)
Vetor solverLU(const Matriz &L, const Matriz &U, const Vetor &b);

// Calcula a inversa de A usando LU
Matriz calcularInversa(const Matriz &A);

#endif