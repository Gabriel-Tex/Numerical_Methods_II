#ifndef PVC_UTILS_H
#define PVC_UTILS_H

#include <vector>

// Resolve sistema tridiagonal Ax = b pelo algoritmo de Thomas
// a: subdiagonal, b: diagonal principal, c: superdiagonal, d: lado direito
std::vector<double> solveTridiagonal(
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> c,
    std::vector<double> d);

// Resolve sistema linear geral Ax = b por eliminação de Gauss com pivotamento
std::vector<double> solveGauss(
    std::vector<std::vector<double>> A,
    std::vector<double> b);

// Solução exata do PVC1: y(x) = 1/(e^{-1} - e) * (e^{-x} - e^x)
double exactSolutionPVC1(double x);

#endif