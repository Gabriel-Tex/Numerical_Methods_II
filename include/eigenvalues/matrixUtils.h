#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <vector>
#include <string>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Retorna a matriz identidade n x n
Matrix identityMatrix(int n);

// Retorna a transposta de A
Matrix transpose(const Matrix &A);

// Produto matriz-vetor: A * v
Vector matVecProduct(const Matrix &A, const Vector &v);

// Produto de duas matrizes: A * B
Matrix matMatProduct(const Matrix &A, const Matrix &B);

// Norma euclidiana de um vetor
double vecNorm(const Vector &v);

// Produto interno de dois vetores
double dotProduct(const Vector &a, const Vector &b);

// Normaliza um vetor (divide pela norma)
Vector normalizeVec(const Vector &v);

// Soma dos quadrados dos elementos abaixo da diagonal principal
double sumSquaresBelowDiagonal(const Matrix &A, int n);

// Imprime uma matriz formatada
void printMatrix(const Matrix &A, const std::string &name = "");

// Imprime um vetor formatado
void printVector(const Vector &v, const std::string &name = "");

// Decomposição LU sem pivotamento: A = L * U
void luDecomposition(const Matrix &A, Matrix &L, Matrix &U);

// Resolve sistema triangular inferior L*y = b
Vector solveLower(const Matrix &L, const Vector &b);

// Resolve sistema triangular superior U*x = y
Vector solveUpper(const Matrix &U, const Vector &y);

// Resolve A*x = b usando decomposição LU pré-calculada
Vector solveLU(const Matrix &L, const Matrix &U, const Vector &b);

// Calcula a inversa de A via LU
Matrix invertMatrix(const Matrix &A);

#endif