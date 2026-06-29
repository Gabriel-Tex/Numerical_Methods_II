#include "../../include/eigenvalues/menuEigenvalues.h"
#include "../../include/eigenvalues/matrixUtils.h"
#include "../../include/eigenvalues/powerMethod.h"
#include "../../include/eigenvalues/householder.h"
#include "../../include/eigenvalues/jacobiMethod.h"
#include "../../include/eigenvalues/qrMethod.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

// ---- utilitários locais ----

static void pause()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int readInt(const std::string &prompt, int min, int max)
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

static double readDouble(const std::string &prompt)
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

static Matrix readMatrix(int n)
{
    Matrix A(n, Vector(n));
    std::cout << "Digite os elementos linha por linha:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            std::cout << "  A[" << i + 1 << "][" << j + 1 << "] = ";
            std::cin >> A[i][j];
        }
    return A;
}

static Vector readVector(int n, const std::string &name)
{
    Vector v(n);
    std::cout << "Digite os elementos de " << name << ":\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << "  v[" << i + 1 << "] = ";
        std::cin >> v[i];
    }
    return v;
}

// Matriz das Tarefas 11-15
static Matrix taskMatrix()
{
    return {{40, 8, 4, 2, 1},
            {8, 30, 12, 6, 2},
            {4, 12, 20, 1, 2},
            {2, 6, 1, 25, 4},
            {1, 2, 2, 4, 5}};
}

// Imprime pares (autovalor, autovetor)
static void printEigenPairs(const Vector &eigenvalues, const Matrix &P)
{
    int n = (int)eigenvalues.size();
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < n; i++)
    {
        std::cout << "\n  Autovalor[" << i + 1 << "] = " << eigenvalues[i] << "\n";
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

// ---- Métodos da Potêcia ----
static void menuPowerMethod()
{
    std::cout << "\n=== Metodos da Potencia ===\n\n";
    std::cout << "Metodo:\n";
    std::cout << "  1 - Potencia Regular   (autovalor dominante)\n";
    std::cout << "  2 - Potencia Inverso   (autovalor de menor |lambda|)\n";
    std::cout << "  3 - Potencia com Deslocamento (autovalor proximo de mu)\n";
    int method = readInt("Escolha: ", 1, 3);

    std::cout << "\nMatriz:\n";
    std::cout << "  1 - Tarefa 11: A1 (3x3)  [[5,2,1],[2,3,1],[1,1,2]]\n";
    std::cout << "  2 - Tarefa 12: A1 (3x3)  [[5,2,1],[2,3,1],[1,1,2]]\n";
    std::cout << "  3 - Tarefa 12: A2 (3x3)  [[-14,1,-2],[1,-1,1],[-2,1,-11]]\n";
    std::cout << "  4 - Tarefa 11/12: A3 (5x5) (matriz das tarefas)\n";
    std::cout << "  5 - Inserir manualmente\n";
    int choice = readInt("Escolha: ", 1, 5);

    Matrix A;
    int n;

    if (choice == 1 || choice == 2)
    {
        n = 3;
        A = {{5, 2, 1}, {2, 3, 1}, {1, 1, 2}};
        std::cout << "\nMatriz A1 (3x3) selecionada.\n";
    }
    else if (choice == 3)
    {
        n = 3;
        A = {{-14, 1, -2}, {1, -1, 1}, {-2, 1, -11}};
        std::cout << "\nMatriz A2 (3x3) selecionada.\n";
    }
    else if (choice == 4)
    {
        n = 5;
        A = taskMatrix();
        std::cout << "\nMatriz A3 (5x5) selecionada.\n";
    }
    else
    {
        n = readInt("Dimensao n da matriz: ", 1, 20);
        A = readMatrix(n);
    }

    printMatrix(A, "A");

    std::cout << "\nVetor inicial v0:\n";
    std::cout << "  1 - Vetor de uns automaticamente\n";
    std::cout << "  2 - Inserir manualmente\n";
    int vChoice = readInt("Escolha: ", 1, 2);

    Vector v0(n, 1.0);
    if (vChoice == 2)
        v0 = readVector(n, "v0");

    double eps = readDouble("Tolerancia eps (ex: 1e-6): ");

    double lambda;
    Vector x;

    if (method == 1)
    {
        regularPower(A, v0, eps, lambda, x);
        std::cout << "\n--- Resultado: Potencia Regular ---\n";
    }
    else if (method == 2)
    {
        inversePower(A, v0, eps, lambda, x);
        std::cout << "\n--- Resultado: Potencia Inverso ---\n";
    }
    else
    {
        double mu = readDouble("Deslocamento mu: ");
        shiftedPower(A, v0, eps, mu, lambda, x);
        std::cout << "\n--- Resultado: Potencia com Deslocamento (mu=" << mu << ") ---\n";
    }

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  Autovalor  lambda = " << lambda << "\n";
    printVector(x, "  Autovetor x");

    pause();
}

// ---- Householder ----
static void menuHouseholder()
{
    std::cout << "\n=== Metodo de Householder ===\n";
    std::cout << "Tarefa 13: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matrix A = taskMatrix();
    printMatrix(A, "A");

    Matrix Abar, H;
    householderMethod(A, Abar, H);

    std::cout << "\n--- Item 1i: Matriz tridiagonal A_bar ---\n";
    printMatrix(Abar, "A_bar");

    std::cout << "\n--- Item 1ii: Matriz acumulada H = H1*H2*H3 ---\n";
    printMatrix(H, "H");

    std::cout << "\n--- Item 3: Autovalores de A_bar via Potencia com Deslocamento ---\n";

    std::vector<double> mus = {40.0, 28.0, 19.0, 5.0, 3.0};
    double eps = 1e-8;
    Vector v0(n, 1.0);

    std::vector<double> eigenvalues;
    std::vector<Vector> eigenvectorsAbar;

    for (int k = 0; k < n; k++)
    {
        double lk;
        Vector xk;
        shiftedPower(Abar, v0, eps, mus[k], lk, xk);
        eigenvalues.push_back(lk);
        eigenvectorsAbar.push_back(xk);
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  Autovalor[" << k + 1 << "] de A_bar = " << lk << "\n";
        printVector(xk, "  Autovetor de A_bar");
    }

    std::cout << "\n--- Item 4: Autovetores de A = H * v_bar ---\n";
    for (int k = 0; k < n; k++)
    {
        Vector xA = matVecProduct(H, eigenvectorsAbar[k]);
        xA = normalizeVec(xA);
        std::cout << "  Autovalor[" << k + 1 << "] = " << eigenvalues[k] << "\n";
        printVector(xA, "  Autovetor de A");
    }

    std::cout << "\n--- Item 5: Autovalores de A ---\n";
    for (int k = 0; k < n; k++)
        std::cout << "  lambda[" << k + 1 << "] = " << eigenvalues[k] << "\n";

    pause();
}

// ---- Jacobi ----
static void menuJacobi()
{
    std::cout << "\n=== Metodo de Jacobi ===\n";
    std::cout << "Tarefa 14: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matrix A = taskMatrix();
    printMatrix(A, "A");

    double eps = 1e-10;

    std::cout << "\n--- Item 1iii: Matrizes apos cada varredura de Jacobi ---\n";

    Matrix Aprev = A;
    Matrix P = identityMatrix(n);
    double val = 100.0;
    int sweep = 0;

    while (val > eps)
    {
        Matrix Anew, Pvar;
        jacobiSweep(Aprev, n, Anew, Pvar);
        Aprev = Anew;
        P = matMatProduct(P, Pvar);
        val = sumSquaresBelowDiagonal(Anew, n);
        sweep++;

        std::cout << "\nVarredura " << sweep << " (soma^2 abaixo diag = "
                  << std::scientific << val << "):\n";
        printMatrix(Anew);
    }

    Vector eigenvalues(n);
    for (int i = 0; i < n; i++)
        eigenvalues[i] = Aprev[i][i];

    std::cout << "\n--- Item 1i: Matriz diagonal A_bar ---\n";
    printMatrix(Aprev, "A_bar (diagonal)");

    std::cout << "\n--- Item 1ii: Matriz acumulada P ---\n";
    printMatrix(P, "P");

    std::cout << "\n--- Item 1iv: Pares (autovalor, autovetor) de A ---\n";
    printEigenPairs(eigenvalues, P);

    std::cout << "\n\n--- Item 2: Jacobi sobre a saida do Householder ---\n";
    Matrix AbarH, H;
    householderMethod(A, AbarH, H);
    std::cout << "Matriz tridiagonal de Householder:\n";
    printMatrix(AbarH);

    std::cout << "\nAplicando Jacobi sobre A_bar (tridiagonal)...\n";
    Matrix P2;
    Vector eigenvalues2;
    jacobiMethod(AbarH, n, eps, P2, eigenvalues2);

    std::cout << "\n--- Item 2ii: As colunas de P2 NAO sao autovetores de A ---\n";
    printMatrix(P2, "P2 (autovetores de A_bar, nao de A)");

    Matrix X = matMatProduct(H, P2);
    std::cout << "\n--- Item 2iv: X = H * P2 (colunas sao autovetores de A) ---\n";
    printEigenPairs(eigenvalues2, X);

    pause();
}

// ---- QR ----
static void menuQR()
{
    std::cout << "\n=== Metodo QR ===\n";
    std::cout << "Tarefa 15: Matriz A (5x5) da tarefa\n\n";

    int n = 5;
    Matrix A = taskMatrix();
    printMatrix(A, "A");

    double eps = 1e-10;

    std::cout << "\n--- Item 1iii: Matrizes apos cada iteracao QR ---\n";

    Matrix Aprev = A;
    Matrix P = identityMatrix(n);
    double val = 100.0;
    int iter = 0;

    while (val > eps)
    {
        Matrix Q, R;
        qrDecomposition(Aprev, n, Q, R);
        Matrix Anew = matMatProduct(R, Q);
        Aprev = Anew;
        P = matMatProduct(P, Q);
        val = sumSquaresBelowDiagonal(Anew, n);
        iter++;

        if (iter <= 3 || val <= eps * 100)
        {
            std::cout << "\nIteracao " << iter << " (soma^2 abaixo diag = "
                      << std::scientific << val << "):\n";
            printMatrix(Anew);
        }
    }
    std::cout << "\nTotal de iteracoes QR: " << iter << "\n";

    Vector eigenvalues(n);
    for (int i = 0; i < n; i++)
        eigenvalues[i] = Aprev[i][i];

    std::cout << "\n--- Item 1i: Matriz diagonal A_bar ---\n";
    printMatrix(Aprev, "A_bar (diagonal)");

    std::cout << "\n--- Item 1ii: Matriz acumulada P ---\n";
    printMatrix(P, "P");

    std::cout << "\n--- Item 1iv: Pares (autovalor, autovetor) de A ---\n";
    printEigenPairs(eigenvalues, P);

    std::cout << "\n\n--- Item 2: QR sobre a saida do Householder ---\n";
    Matrix AbarH, H;
    householderMethod(A, AbarH, H);
    std::cout << "Matriz tridiagonal de Householder:\n";
    printMatrix(AbarH);

    std::cout << "\nAplicando QR sobre A_bar (tridiagonal)...\n";

    Matrix A2prev = AbarH;
    Matrix P2 = identityMatrix(n);
    val = 100.0;
    iter = 0;

    while (val > eps)
    {
        Matrix Q, R;
        qrDecomposition(A2prev, n, Q, R);
        Matrix A2new = matMatProduct(R, Q);
        A2prev = A2new;
        P2 = matMatProduct(P2, Q);
        val = sumSquaresBelowDiagonal(A2new, n);
        iter++;

        if (iter == 1)
        {
            std::cout << "\n--- Item 2i: Apos 1a iteracao QR sobre tridiagonal ---\n";
            printMatrix(A2new);
        }
    }
    std::cout << "\nTotal de iteracoes QR sobre tridiagonal: " << iter << "\n";

    Vector eigenvalues2(n);
    for (int i = 0; i < n; i++)
        eigenvalues2[i] = A2prev[i][i];

    std::cout << "\n--- Item 2ii: Colunas de P2 NAO sao autovetores de A ---\n";
    printMatrix(P2, "P2 (autovetores de A_bar, nao de A)");

    Matrix X = matMatProduct(H, P2);
    std::cout << "\n--- Item 2iv: X = H * P2 (colunas sao autovetores de A) ---\n";
    printEigenPairs(eigenvalues2, X);

    pause();
}

// ---- Menu principal de autovalores ----
void menuEigenvalues()
{
    while (true)
    {
        std::cout << "   Autovalores e Autovetores              \n";
        std::cout << "  1 - Metodos da Potencia                 \n";
        std::cout << "  2 - Householder                         \n";
        std::cout << "  3 - Jacobi                              \n";
        std::cout << "  4 - QR                                  \n";
        std::cout << "  0 - Voltar                              \n";

        int op = readInt("Opcao: ", 0, 4);
        if (op == 0)
            return;
        if (op == 1)
            menuPowerMethod();
        if (op == 2)
            menuHouseholder();
        if (op == 3)
            menuJacobi();
        if (op == 4)
            menuQR();
    }
}