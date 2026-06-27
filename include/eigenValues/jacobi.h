#ifndef JACOBI_H
#define JACOBI_H

#include "matrizUtils.h"

//  Método de Jacobi
//
// Diagonaliza a matriz simétrica A através de múltiplas
// varreduras de Jacobi, até que os elementos fora da diagonal
// sejam menores que eps.
//
// Retorna:
//   P - matriz acumulada dos produtos das matrizes de Jacobi (suas colunas são os autovetores de A)
//   Lamb - vetor com os autovalores de A (diagonal de A_barra)
void metodoDeJacobi(const Matriz &A,
                    int n,
                    double eps,
                    Matriz &P,
                    Vetor &Lamb);

// Executa uma varredura de Jacobi sobre a matriz A.
// Percorre todas as colunas j=0..n-2 e linhas i=j+1..n-1,
// zerando os elementos fora da diagonal via transformações
// de similaridade com matrizes de rotação de Jacobi.
//
// Retorna:
//   A_nova - matriz após a varredura
//   P_var - produto das matrizes de Jacobi desta varredura
void varreduraDeJacobi(const Matriz &A,
                       int n,
                       Matriz &A_nova,
                       Matriz &P_var);

// Constrói a matriz de rotação de Jacobi G_ij para zerar
// o elemento (i, j) de A via transformação de similaridade.
// Baseado na seção 3.1.4.2 das notas de aula.
// i e j são 0-based.
Matriz matrizJacobiBaseadaNoElemento_ij_DaMatrizVelha(const Matriz &A,
                                                      int i, int j, int n);

#endif