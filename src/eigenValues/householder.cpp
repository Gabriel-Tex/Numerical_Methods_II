#include "../../include/eigenValues/householder.h"
#include <cmath>

// Construção da matriz de Householder H_i
//
// Dada a matriz A_{i-1} e o índice i (1-based conforme notas),
// constrói H_i que zerará os elementos abaixo de (i+1, i)
// na coluna i, preservando a estrutura tridiagonal já obtida.
Matriz matrizHouseholderBaseadaNaCol(const Matriz &A, int i)
{
    // i é 1-based
    // internamente usamos índices 0-based
    int n = (int)A.size();

    // w: vetor nulo com n elementos
    // Copia os elementos abaixo da diagonal da coluna i
    // para as posições (i+1) até (n-1) - índices 0-based: i até n-1
    Vetor w(n, 0.0);
    for (int k = i; k < n; k++) // posição i-1+1 = i em 0-based
        w[k] = A[k][i - 1]; // coluna i em 1-based = i-1 em 0-based

    // Calcular o comprimento de w
    double Lw = normaVetor(w);

    // w': copia Lw na posição i (0-based)
    Vetor wPrime(n, 0.0);
    wPrime[i] = Lw; // posição i+1 em 1-based = i em 0-based

    // N = w - w'
    Vetor N(n, 0.0);
    for (int k = 0; k < n; k++)
        N[k] = w[k] - wPrime[k];

    // ||N||
    double LN = normaVetor(N);

    // Se N é nulo, retorna identidade (coluna já zerada)
    if (LN < 1e-14)
        return matrizIdentidade(n);

    // n_vec = N / ||N||
    Vetor n_vec(n);
    for (int k = 0; k < n; k++)
        n_vec[k] = N[k] / LN;

    // Montar H = I - 2 * n_vec * n_vec^T
    Matriz Hi = matrizIdentidade(n);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            Hi[r][c] -= 2.0 * n_vec[r] * n_vec[c];

    return Hi;
}

//  Método de Householder 
//
// Loop de i = 1 até n-2 (1-based), ou seja n-2 passos.
void metodoDeHouseholder(const Matriz &A,
                         Matriz &A_barra,
                         Matriz &H)
{
    int n = (int)A.size();

    // H ← I  (acumulará o produto I * H1 * H2 * ... * H_{n-2})
    H = matrizIdentidade(n);

    // A_{i-1} ← A
    Matriz A_prev = A;

    // Loop i = 1 ... (n-2)  
    for (int i = 1; i <= n - 2; i++)
    {
        // Construir H_i baseada na coluna i de A_{i-1}
        Matriz Hi = matrizHouseholderBaseadaNaCol(A_prev, i);

        // Transformação de similaridade: A_i = H_i^T * A_{i-1} * H_i
        // Como H_i é simétrica (H_i^T = H_i), fica: H_i * A_{i-1} * H_i
        Matriz Ai = produtoMatrizes(produtoMatrizes(Hi, A_prev), Hi);

        // Salvar para o próximo passo
        A_prev = Ai;

        // Acumular: H ← H * H_i
        H = produtoMatrizes(H, Hi);
    }

    A_barra = A_prev;
}