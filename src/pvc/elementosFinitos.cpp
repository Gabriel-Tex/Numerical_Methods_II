#include "../../include/pvc/elementosFinitos.h"
#include "../../include/pvc/pvcUtils.h"
#include <cmath>

// ------------------------------------------------------------------
// MEF para PVC1 (Aula 28, secao 4.1)
//
// Formulacao fraca apos integracao por partes (eq. 16):
//   integral(y'^w') + integral(y*w) = 0
//
// Funcoes de interpolacao lineares em cada elemento.
// Matriz local 2x2 de cada elemento (eq. 24):
//   K_e = (1/h)*[ 1 -1; -1  1] + (h/6)*[2 1; 1 2]
//
// Nos interiores: 1 .. N-1 (M = N-1 incognitas)
// ------------------------------------------------------------------
std::vector<double> elementosFinitasPVC1(int N)
{
    int M = N - 1; // nos interiores (incognitas)
    double h = 1.0 / N;

    // Contribuicao do elemento (eq. 24):
    //   rigidez: (1/h)*[ 1 -1; -1  1]
    //   massa:   (h/6)*[ 2  1;  1  2]
    double K11 = 1.0/h + 2.0*h/6.0; // diagonal da matriz local
    double K12 = -1.0/h + h/6.0;    // fora da diagonal da matriz local

    // Montagem: matriz global tridiagonal de tamanho M x M
    // Nos interiores sao indexados 1..N-1, mapeados para 0..M-1
    std::vector<double> a(M, K12);  // subdiagonal
    std::vector<double> b(M, 0.0); // diagonal
    std::vector<double> c(M, K12); // superdiagonal
    std::vector<double> d(M, 0.0);

    // Cada no interior i recebe contribuicao de dois elementos
    for (int k = 0; k < M; k++)
        b[k] = 2.0 * K11; // elemento da esquerda + elemento da direita

    // Nos dos extremos (0 e M-1 do sistema, ie nos 1 e N-1 do dominio)
    // recebem contribuicao de apenas um elemento extra do contorno.
    // Mas o no 0 (x=0) tem Dirichlet y=0, e o no N (x=1) tem y=1.
    // Correcao: nos interiores 1 e N-1 estao no limite, entao
    // sua diagonal recebe apenas dois K11 (de dois elementos adjacentes).
    // Isso ja esta correto acima. O que muda e o lado direito:
    // O elemento mais proximo do no N contribui com -K12 * y_N = -K12 * 1.
    d[M - 1] = -K12 * 1.0; // contribuicao de y_{N} = 1 (Dirichlet)

    // Correcao do primeiro e ultimo no: recebem so um elemento
    // No 1 (indice 0): elemento el_1 contribui com K11 do lado direito (no de contorno y_0=0)
    // -> contribuicao e zero (y_0 = 0)
    // No N-1 (indice M-1): elemento el_N contribui com K11 a diagonal
    // ja esta somado. Para nos nos extremos do dominio interior,
    // a diagonal recebe K11 apenas de um elemento.
    // Corrigindo: o primeiro e ultimo no interior recebem K11 de apenas UM elemento
    // do lado do contorno, mais K11 do elemento vizinho interior.
    // A formula 2*K11 ja esta correta porque cada no interior esta entre dois elementos.
    // (o no 1 esta entre el_1 e el_2, o no N-1 esta entre el_{N-1} e el_N)

    return resolveTridiagonal(a, b, c, d);
}

// ------------------------------------------------------------------
// MEF para PVC2 (Aula 28, secao 4.2)
//
// Elemento quadrilateral bilinear 2D, 4 nos.
// Funcoes de interpolacao de Legendre no elemento padrao [-1,1]^2 (eq 30-33).
//
// Matriz local 4x4 de coeficientes (eq. 55, resultado da Aula 28):
//   Para elementos quadrados dx=dy=1/N, a matriz e:
//
//   Ke = [ 2/3  -1/6  -1/3  -1/6 ]
//        [-1/6   2/3  -1/6  -1/3 ]
//        [-1/3  -1/6   2/3  -1/6 ]
//        [-1/6  -1/3  -1/6   2/3 ]
//
// Vetor local de carga (eq. 56, com f=4):
//   fe = -1/16 (cada componente, pois f*|J|*integral = 4*(1/64)*4 = 1/4... ver calculo)
//
// Numeracao dos nos do elemento (Aula 28, convencao):
//   no local 1: canto inferior esquerdo (i,   j  )
//   no local 2: canto inferior direito  (i+1, j  )
//   no local 3: canto superior direito  (i+1, j+1)
//   no local 4: canto superior esquerdo (i,   j+1)
// ------------------------------------------------------------------
std::vector<double> elementosFinitasPVC2(int N)
{
    int M = N - 1;   // nos interiores por direcao
    int total = M * M;

    // Matriz local 4x4 (Aula 28, eq. 55)
    // Calculada analiticamente para elementos quadrados de lado h=1/N
    double a11 =  2.0/3.0;
    double a12 = -1.0/6.0;
    double a13 = -1.0/3.0;
    double a14 = -1.0/6.0;
    // Por simetria:
    // a21=a12, a22=a11, a23=a14, a24=a13
    // a31=a13, a32=a14, a33=a11, a34=a12
    // a41=a14, a42=a13, a43=a12, a44=a11
    double Ke[4][4] = {
        { a11, a12, a13, a14 },
        { a12, a11, a14, a13 },
        { a13, a14, a11, a12 },
        { a14, a13, a12, a11 }
    };

    // Vetor de carga local (Aula 28, eq. 56): fe_i = -1/16
    double fe = -1.0/16.0;

    // Sistema global
    std::vector<std::vector<double>> A(total, std::vector<double>(total, 0.0));
    std::vector<double> rhs(total, 0.0);

    // Numeracao dos nos globais: no (i,j) -> indice j*(N+1) + i
    // Nos interiores: i=1..N-1, j=1..N-1
    // Mapeamento para indice do sistema: (i-1) + (j-1)*M
    auto gidx = [&](int i, int j) -> int {
        // i, j sao indices globais de nos (0..N)
        // Retorna -1 se no do contorno, ou indice do sistema (0..M^2-1)
        if (i == 0 || i == N || j == 0 || j == N) return -1;
        return (j - 1) * M + (i - 1);
    };

    // Itera sobre todos os elementos
    for (int ej = 0; ej < N; ej++)       // indice do elemento em y
    {
        for (int ei = 0; ei < N; ei++)   // indice do elemento em x
        {
            // Nos globais do elemento (convencao Aula 28):
            // no local 1: (ei,   ej  )
            // no local 2: (ei+1, ej  )
            // no local 3: (ei+1, ej+1)
            // no local 4: (ei,   ej+1)
            int gn[4] = {
                gidx(ei,   ej  ),
                gidx(ei+1, ej  ),
                gidx(ei+1, ej+1),
                gidx(ei,   ej+1)
            };

            // Montagem
            for (int a = 0; a < 4; a++)
            {
                if (gn[a] < 0) continue; // no de contorno: Dirichlet=0, sem contribuicao
                rhs[gn[a]] += fe;        // contribuicao do vetor de carga

                for (int b = 0; b < 4; b++)
                {
                    if (gn[b] < 0) continue; // no de contorno (Dirichlet=0)
                    A[gn[a]][gn[b]] += Ke[a][b];
                }
            }
        }
    }

    return resolveGauss(A, rhs);
}