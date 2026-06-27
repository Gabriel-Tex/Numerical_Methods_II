#include "../../include/eigenValues/potencia.h"
#include <cmath>

// ========================== Método da Potência Regular ==========================
void potenciaRegular(const Matriz &A,
                     const Vetor &v0,
                     double eps,
                     double &lambda,
                     Vetor &x)
{
    int n = (int)A.size();

    // Step 2: inicializar lambda
    double lambdaNovo = 0.0;

    // Step 3: copiar v0 para vk_Novo
    Vetor vk_Novo = v0;

    double lambdaVelho;
    Vetor x1_Velho(n);

    do
    {
        // Step 4: copiar lambdaNovo para lambdaVelho
        lambdaVelho = lambdaNovo;

        // Step 5: copiar vk_Novo para vk_Velho
        Vetor vk_Velho = vk_Novo;

        // Step 6: normalizar vk_Velho → x1_Velho
        x1_Velho = normalizaVetor(vk_Velho);

        // Step 7: calcular vk_Novo = A * x1_Velho
        vk_Novo = produtoMatrizVetor(A, x1_Velho);

        // Step 8: nova estimativa de lambda = x1_Velho^T * vk_Novo
        lambdaNovo = produtoInterno(x1_Velho, vk_Novo);

        // Step 9: verificar convergência
    } while (std::fabs((lambdaNovo - lambdaVelho) / lambdaNovo) > eps);

    lambda = lambdaNovo;
    x = x1_Velho;
}

// ========================== Aula 19 — Método da Potência Inverso ==========================
void potenciaInverso(const Matriz &A,
                     const Vetor &v0,
                     double eps,
                     double &lambda,
                     Vetor &x)
{
    int n = (int)A.size();

    // Step 2: calcular decomposição LU de A
    Matriz L, U;
    decompLU(A, L, U);

    // Step 3: inicializar lambda_bar
    double lambdaBarNovo = 0.0;

    // Step 4: copiar v0 para vk_Novo
    Vetor vk_Novo = v0;

    double lambdaBarVelho;
    Vetor x1_Velho(n);

    do
    {
        // Step 5
        lambdaBarVelho = lambdaBarNovo;

        // Step 6: copiar vk_Novo para vk_Velho
        Vetor vk_Velho = vk_Novo;

        // Step 7: normalizar vk_Velho → x1_Velho
        x1_Velho = normalizaVetor(vk_Velho);

        // Step 8: calcular vk_Novo = A^{-1} * x1_Velho usando LU
        vk_Novo = solverLU(L, U, x1_Velho);

        // Step 9: nova estimativa de lambda_bar = x1_Velho^T * vk_Novo
        lambdaBarNovo = produtoInterno(x1_Velho, vk_Novo);

        // Step 10: verificar convergência
    } while (std::fabs((lambdaBarNovo - lambdaBarVelho) / lambdaBarNovo) > eps);

    // Step 11: lambda_n = 1 / lambda_bar
    lambda = 1.0 / lambdaBarNovo;

    // Step 12: autovetor é x1_Velho
    x = x1_Velho;
}

// ========================== Método da Potência com Deslocamento ==========================
void potenciaComDeslocamento(const Matriz &A,
                             const Vetor &v0,
                             double eps,
                             double mu,
                             double &lambda,
                             Vetor &x)
{
    int n = (int)A.size();

    // Step 1: A_hat = A - mu * I
    Matriz A_hat = A;
    for (int i = 0; i < n; i++)
        A_hat[i][i] -= mu;

    // Step 2: aplicar potência inverso sobre A_hat
    double lambdaHat;
    Vetor xHat;
    potenciaInverso(A_hat, v0, eps, lambdaHat, xHat);

    // Step 3: lambda_i = lambda_hat + mu
    lambda = lambdaHat + mu;

    // Step 4: autovetor é o mesmo
    x = xHat;
}