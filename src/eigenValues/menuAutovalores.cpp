#include "../../include/eigenValues/menuAutovalores.h"
#include "../../include/eigenValues/matrizUtils.h"
#include "../../include/eigenValues/potencia.h"
#include "../../include/eigenValues/householder.h"
#include "../../include/eigenValues/jacobi.h"
#include "../../include/eigenValues/metodoQR.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

// ---- Utilitários locais ----

static void pausarEigen()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int lerIntEigen(const std::string &prompt, int min, int max)
{
    int v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v && v >= min && v <= max)
            return v;
        std::cout << "  Opcao invalida. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static double lerDoubleEigen(const std::string &prompt)
{
    double v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v)
            return v;
        std::cout << "  Valor invalido. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static Matriz lerMatriz(int n)
{
    Matriz A(n, Vetor(n));
    std::cout << "Digite os elementos linha por linha:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            std::cout << "  A[" << i + 1 << "][" << j + 1 << "] = ";
            std::cin >> A[i][j];
        }
    return A;
}

static Vetor lerVetor(int n, const std::string &nome)
{
    Vetor v(n);
    std::cout << "Digite os elementos de " << nome << ":\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << "  v[" << i + 1 << "] = ";
        std::cin >> v[i];
    }
    return v;
}

// ---- Matriz A da Tarefa 11, 12, 13, 14, 15 ----
static Matriz matrizTarefa()
{
    return {{40, 8, 4, 2, 1},
            {8, 30, 12, 6, 2},
            {4, 12, 20, 1, 2},
            {2, 6, 1, 25, 4},
            {1, 2, 2, 4, 5}};
}

// ---- Imprime pares autovalor, autovetor ----
static void imprimirPares(const Vetor &Lamb, const Matriz &P)
{
    int n = (int)Lamb.size();
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < n; i++)
    {
        std::cout << "\n  Autovalor[" << i + 1 << "] = " << Lamb[i] << "\n";
        std::cout << "  Autovetor[" << i + 1 << "] = [";
        for (int k = 0; k < n; k++)
        {
            std::cout << P[k][i];
            if (k + 1 < n)
                std::cout << "  ";
        }
        std::cout << "]\n";
    }
}

//  Métodos da Potência
void menuPotencia()
{
    std::cout << "\n=== Metodos da Potencia ===\n\n";
    std::cout << "Metodo:\n";
    std::cout << "  1 - Potencia Regular   (autovalor dominante)\n";
    std::cout << "  2 - Potencia Inverso   (autovalor de menor |lambda|)\n";
    std::cout << "  3 - Potencia com Deslocamento (autovalor proximo de mu)\n";
    int metodo = lerIntEigen("Escolha: ", 1, 3);

    std::cout << "\nMatriz:\n";
    std::cout << "  1 - Tarefa 11: A1 (3x3)  [[5,2,1],[2,3,1],[1,1,2]]\n";
    std::cout << "  2 - Tarefa 12: A1 (3x3)  [[5,2,1],[2,3,1],[1,1,2]]\n";
    std::cout << "  3 - Tarefa 12: A2 (3x3)  [[-14,1,-2],[1,-1,1],[-2,1,-11]]\n";
    std::cout << "  4 - Tarefa 11/12: A3 (5x5) (matriz das tarefas)\n";
    std::cout << "  5 - Inserir manualmente\n";
    int escolha = lerIntEigen("Escolha: ", 1, 5);

    Matriz A;
    int n;

    if (escolha == 1 || escolha == 2)
    {
        n = 3;
        A = {{5, 2, 1}, {2, 3, 1}, {1, 1, 2}};
        std::cout << "\nMatriz A1 (3x3) selecionada.\n";
    }
    else if (escolha == 3)
    {
        n = 3;
        A = {{-14, 1, -2}, {1, -1, 1}, {-2, 1, -11}};
        std::cout << "\nMatriz A2 (3x3) selecionada.\n";
    }
    else if (escolha == 4)
    {
        n = 5;
        A = matrizTarefa();
        std::cout << "\nMatriz A3 (5x5) selecionada.\n";
    }
    else
    {
        n = lerIntEigen("Dimensao n da matriz: ", 1, 20);
        A = lerMatriz(n);
    }

    imprimirMatriz(A, "A");

    std::cout << "\nVetor inicial v0 (pode usar [1,1,...,1]):\n";
    std::cout << "  1 - Usar vetor de uns automaticamente\n";
    std::cout << "  2 - Inserir manualmente\n";
    int vEscolha = lerIntEigen("Escolha: ", 1, 2);

    Vetor v0(n, 1.0);
    if (vEscolha == 2)
        v0 = lerVetor(n, "v0");

    double eps = lerDoubleEigen("Tolerancia eps (ex: 1e-6): ");

    double lambda;
    Vetor x;

    if (metodo == 1)
    {
        potenciaRegular(A, v0, eps, lambda, x);
        std::cout << "\n--- Resultado: Potencia Regular ---\n";
    }
    else if (metodo == 2)
    {
        potenciaInverso(A, v0, eps, lambda, x);
        std::cout << "\n--- Resultado: Potencia Inverso ---\n";
    }
    else
    {
        double mu = lerDoubleEigen("Deslocamento mu: ");
        potenciaComDeslocamento(A, v0, eps, mu, lambda, x);
        std::cout << "\n--- Resultado: Potencia com Deslocamento (mu=" << mu << ") ---\n";
    }

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  Autovalor  lambda = " << lambda << "\n";
    imprimirVetor(x, "  Autovetor x");

    pausarEigen();
}

//  Método de Householder
void menuHouseholder()
{
    std::cout << "\n=== Metodo de Householder ===\n";
    std::cout << "Tarefa 13: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matriz A = matrizTarefa();
    imprimirMatriz(A, "A");

    // Item 1: aplicar Householder
    Matriz A_barra, H;
    metodoDeHouseholder(A, A_barra, H);

    std::cout << "\n--- Item 1i: Matriz tridiagonal A_barra ---\n";
    imprimirMatriz(A_barra, "A_barra");

    std::cout << "\n--- Item 1ii: Matriz acumulada H = H1*H2*H3 ---\n";
    imprimirMatriz(H, "H");

    // Item 3: usar potencia regular para achar autovalores/vetores de A_barra
    std::cout << "\n--- Item 3: Autovalores de A_barra via Potencia com Deslocamento ---\n";
    std::cout << "(Usando deslocamentos estimados para cada autovalor)\n";

    // Estimativas dos deslocamentos baseadas na diagonal de A_barra
    std::vector<double> mus = {40.0, 28.0, 19.0, 5.0, 3.0};
    double eps = 1e-8;
    Vetor v0(n, 1.0);

    std::vector<double> lambdas;
    std::vector<Vetor> autovetoresAbarra;

    for (int k = 0; k < n; k++)
    {
        double lk;
        Vetor xk;
        potenciaComDeslocamento(A_barra, v0, eps, mus[k], lk, xk);
        lambdas.push_back(lk);
        autovetoresAbarra.push_back(xk);
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  Autovalor[" << k + 1 << "] de A_barra = " << lk << "\n";
        imprimirVetor(xk, "  Autovetor de A_barra");
    }

    // Item 4: autovetores de A = H * autovetores de A_barra
    std::cout << "\n--- Item 4: Autovetores de A = H * v_barra ---\n";
    for (int k = 0; k < n; k++)
    {
        Vetor xA = produtoMatrizVetor(H, autovetoresAbarra[k]);
        xA = normalizaVetor(xA);
        std::cout << "  Autovalor[" << k + 1 << "] = " << lambdas[k] << "\n";
        imprimirVetor(xA, "  Autovetor de A");
    }

    // Item 5: autovalores de A = autovalores de A_barra
    std::cout << "\n--- Item 5: Autovalores de A ---\n";
    for (int k = 0; k < n; k++)
        std::cout << "  lambda[" << k + 1 << "] = " << lambdas[k] << "\n";

    pausarEigen();
}

// Método de Jacobi
void menuJacobi()
{
    std::cout << "\n=== Metodo de Jacobi ===\n";
    std::cout << "Tarefa 14: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matriz A = matrizTarefa();
    imprimirMatriz(A, "A");

    double eps = 1e-10;

    // Item 1: aplicar Jacobi
    Matriz P;
    Vetor Lamb;

    // Rodar com impressão de cada varredura
    std::cout << "\n--- Item 1iii: Matrizes apos cada varredura de Jacobi ---\n";

    Matriz A_prev = A;
    P = matrizIdentidade(n);
    double val = 100.0;
    int varredura = 0;

    while (val > eps)
    {
        Matriz A_nova, P_var;
        varreduraDeJacobi(A_prev, n, A_nova, P_var);
        A_prev = A_nova;
        P = produtoMatrizes(P, P_var);
        val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A_nova, n);
        varredura++;

        std::cout << "\nVarredura " << varredura << " (soma^2 abaixo diag = "
                  << std::scientific << val << "):\n";
        imprimirMatriz(A_nova);
    }

    // Extrair autovalores
    Lamb.resize(n);
    for (int i = 0; i < n; i++)
        Lamb[i] = A_prev[i][i];

    std::cout << "\n--- Item 1i: Matriz diagonal A_barra ---\n";
    imprimirMatriz(A_prev, "A_barra (diagonal)");

    std::cout << "\n--- Item 1ii: Matriz acumulada P ---\n";
    imprimirMatriz(P, "P");

    std::cout << "\n--- Item 1iv: Pares (autovalor, autovetor) de A ---\n";
    imprimirPares(Lamb, P);

    // Item 2: Jacobi sobre a saida do Householder
    std::cout << "\n\n--- Item 2: Jacobi sobre a saida do Householder ---\n";
    Matriz A_barra_H, H;
    metodoDeHouseholder(A, A_barra_H, H);
    std::cout << "Matriz tridiagonal de Householder:\n";
    imprimirMatriz(A_barra_H);

    std::cout << "\nAplicando Jacobi sobre A_barra (tridiagonal)...\n";
    Matriz P2;
    Vetor Lamb2;
    metodoDeJacobi(A_barra_H, n, eps, P2, Lamb2);

    std::cout << "\n--- Item 2ii: As colunas de P2 NAO sao autovetores de A ---\n";
    imprimirMatriz(P2, "P2 (autovetores de A_barra, nao de A)");

    // Item 2iv: X = H * P2  → autovetores de A
    Matriz X = produtoMatrizes(H, P2);
    std::cout << "\n--- Item 2iv: X = H * P2 (colunas sao autovetores de A) ---\n";
    imprimirPares(Lamb2, X);

    pausarEigen();
}

//  Método QR
void menuQR()
{
    std::cout << "\n=== Metodo QR ===\n";
    std::cout << "Tarefa 15: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matriz A = matrizTarefa();
    imprimirMatriz(A, "A");

    double eps = 1e-10;

    // Item 1: aplicar QR com impressão de cada iteração
    std::cout << "\n--- Item 1iii: Matrizes apos cada iteracao QR ---\n";

    Matriz A_prev = A;
    Matriz P = matrizIdentidade(n);
    double val = 100.0;
    int iter = 0;

    while (val > eps)
    {
        Matriz Q, R;
        decomposicaoQR(A_prev, n, Q, R);
        Matriz A_nova = produtoMatrizes(R, Q);
        A_prev = A_nova;
        P = produtoMatrizes(P, Q);
        val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A_nova, n);
        iter++;

        // Imprimir apenas as primeiras e últimas iterações para não poluir
        if (iter <= 3 || val <= eps * 100)
        {
            std::cout << "\nIteracao " << iter << " (soma^2 abaixo diag = "
                      << std::scientific << val << "):\n";
            imprimirMatriz(A_nova);
        }
    }
    std::cout << "\nTotal de iteracoes QR: " << iter << "\n";

    Vetor Lamb(n);
    for (int i = 0; i < n; i++)
        Lamb[i] = A_prev[i][i];

    std::cout << "\n--- Item 1i: Matriz diagonal A_barra ---\n";
    imprimirMatriz(A_prev, "A_barra (diagonal)");

    std::cout << "\n--- Item 1ii: Matriz acumulada P ---\n";
    imprimirMatriz(P, "P");

    std::cout << "\n--- Item 1iv: Pares (autovalor, autovetor) de A ---\n";
    imprimirPares(Lamb, P);

    // Item 2: QR sobre a saida do Householder
    std::cout << "\n\n--- Item 2: QR sobre a saida do Householder ---\n";
    Matriz A_barra_H, H;
    metodoDeHouseholder(A, A_barra_H, H);
    std::cout << "Matriz tridiagonal de Householder:\n";
    imprimirMatriz(A_barra_H);

    std::cout << "\nAplicando QR sobre A_barra (tridiagonal)...\n";

    Matriz A2_prev = A_barra_H;
    Matriz P2 = matrizIdentidade(n);
    val = 100.0;
    iter = 0;

    while (val > eps)
    {
        Matriz Q, R;
        decomposicaoQR(A2_prev, n, Q, R);
        Matriz A2_nova = produtoMatrizes(R, Q);
        A2_prev = A2_nova;
        P2 = produtoMatrizes(P2, Q);
        val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A2_nova, n);
        iter++;

        // Verificar se termos que eram zero deixaram de ser zero
        if (iter == 1)
        {
            std::cout << "\n--- Item 2i: Apos 1a iteracao QR sobre tridiagonal ---\n";
            std::cout << "(observar se termos zerados voltaram a ser nao-zero)\n";
            imprimirMatriz(A2_nova);
        }
    }
    std::cout << "\nTotal de iteracoes QR sobre tridiagonal: " << iter << "\n";

    Vetor Lamb2(n);
    for (int i = 0; i < n; i++)
        Lamb2[i] = A2_prev[i][i];

    std::cout << "\n--- Item 2ii: Colunas de P2 NAO sao autovetores de A ---\n";
    imprimirMatriz(P2, "P2 (autovetores de A_barra, nao de A)");

    // Item 2iv: X = H * P2 → autovetores de A
    Matriz X = produtoMatrizes(H, P2);
    std::cout << "\n--- Item 2iv: X = H * P2 (colunas sao autovetores de A) ---\n";
    imprimirPares(Lamb2, X);

    pausarEigen();
}