#include "../../include/eigenValues/jacobi.h"
#include <cmath>

// Construção da matriz de rotação de Jacobi G_ij
//
// Índices i e j são 0-based.
//
// Zera o elemento (i, j) de A via transformação de similaridade
// G_ij^T * A_velha * G_ij.
Matriz matrizJacobiBaseadaNoElemento_ij_DaMatrizVelha(const Matriz &A,
                                                      int i, int j, int n)
{
    const double eps = 1e-9;

    // Inicializar Gij com a identidade
    Matriz Gij = matrizIdentidade(n);

    // Se o elemento já é praticamente zero, retorna identidade
    if (std::fabs(A[i][j]) <= eps)
        return Gij;

    double theta;

    // Calcular theta 
    if (std::fabs(A[i][i] - A[j][j]) <= eps)
    {
        // A[i][i] == A[j][j]  → θ = π/4
        theta = M_PI / 4.0;
    }
    else
    {
        // tg(2θ) = -2*A[i][j] / (A[i][i] - A[j][j])
        // θ = (1/2) * arctan(...)
        theta = 0.5 * std::atan(-2.0 * A[i][j] / (A[i][i] - A[j][j]));
    }

    double c = std::cos(theta);
    double s = std::sin(theta);

    // Posições da matriz de rotação
    // el[j][j] = el[i][i] = cos(θ)
    // el[i][j] = sen(θ)
    // el[j][i] = -sen(θ)
    Gij[j][j] = c;
    Gij[i][i] = c;
    Gij[i][j] = s;
    Gij[j][i] = -s;

    return Gij;
}

//  Varredura de Jacobi
void varreduraDeJacobi(const Matriz &A,
                       int n,
                       Matriz &A_nova,
                       Matriz &P_var)
{
    // P_var acumula o produto das matrizes de Jacobi desta varredura
    P_var = matrizIdentidade(n);

    Matriz A_prev = A;

    // Loop j = 0 ... n-2 (colunas)
    for (int j = 0; j <= n - 2; j++)
    {
        // Loop i = j+1 ... n-1 (linhas abaixo da diagonal)
        for (int i = j + 1; i <= n - 1; i++)
        {
            // Construir Gij para zerar A_prev[i][j]
            Matriz Gij = matrizJacobiBaseadaNoElemento_ij_DaMatrizVelha(A_prev, i, j, n);

            // Transformação de similaridade: A_nova = Gij^T * A_prev * Gij
            Matriz GijT = transposta(Gij);
            A_nova = produtoMatrizes(produtoMatrizes(GijT, A_prev), Gij);

            // Salvar para o próximo passo
            A_prev = A_nova;

            // Acumular: P_var ← P_var * Gij
            P_var = produtoMatrizes(P_var, Gij);
        }
    }

    A_nova = A_prev;
}

//  Método de Jacobi (loop de varreduras)
void metodoDeJacobi(const Matriz &A,
                    int n,
                    double eps,
                    Matriz &P,
                    Vetor &Lamb)
{
    // P ← I  (acumulará todas as matrizes de Jacobi)
    P = matrizIdentidade(n);

    Matriz A_prev = A;
    double val = 100.0; // valor inicial para entrar no loop

    while (val > eps)
    {
        Matriz A_nova, P_var;
        varreduraDeJacobi(A_prev, n, A_nova, P_var);

        // Salvar para a próxima varredura
        A_prev = A_nova;

        // Acumular: P ← P * P_var
        P = produtoMatrizes(P, P_var);

        // Verificar convergência
        val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A_nova, n);
    }

    // Extrair autovalores da diagonal de A_prev
    Lamb.resize(n);
    for (int i = 0; i < n; i++)
        Lamb[i] = A_prev[i][i];
}