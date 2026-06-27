#ifndef MULTIDIMINTEGRAL_H
#define MULTIDIMINTEGRAL_H

#include <functional>

using Func2D = std::function<double(double, double)>;

// INTEGRAIS MULTIDIMENSIONAIS
//
// Área de superfície 3D: A = ∫∫ √((∂f/∂x)²+(∂f/∂y)²+1) dA
// Volume abaixo de superfície: V = ∫∫ f(x,y) dA
//
// Tarefa Aula 15 - Problema 2 (domínio elíptico):
//   elipse x²/1600 + y²/1600 ≤ 1 (a=b=40 para superfície)
//   Mudança 1: coordenadas polares elípticas (α∈[0,1], β∈[0,2π])
//   Mudança 2: GL mapeia [0,1]→[-1,1] e [0,2π]→[-1,1]
//
// Tarefa Aula 16 - Problema 2 (domínio elíptico):
//   elipse x²/1600 + y²/400 ≤ 1 (a=40, b=20 para volume)

// Área de superfície 3D sobre domínio retangular [xa,xb]×[ya,yb]
// Usa GL 2D com n pontos em cada direção
double areaSuperficieRetangular(Func2D dfx, Func2D dfy,
                                double xa, double xb, double ya, double yb, int n);

// Área de superfície 3D sobre domínio elíptico (semieixos a, b)
// via mudança de variável polar elíptica + GL com n pontos
double areaSuperficieEliptica(Func2D dfx, Func2D dfy, double a, double b, int n);

// Volume abaixo de superfície sobre domínio retangular
double volumeRetangular(Func2D f, double xa, double xb, double ya, double yb, int n);

// Volume abaixo de superfície sobre domínio elíptico (semieixos a, b)
double volumeEliptico(Func2D f, double a, double b, int n);

#endif