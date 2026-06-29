#include "../../include/pvc/menuPVC.h"
#include "../../include/pvc/finiteDiff.h"
#include "../../include/pvc/finiteElements.h"
#include "../../include/pvc/pvcUtils.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

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

// ---- exibição de tabelas ----

static void printPVC1Table(const std::string &methodName, int N,
                           const std::vector<double> &approx)
{
    int M = N - 1;
    double h = 1.0 / N;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n--- " << methodName << " (N=" << N << ") ---\n";
    std::cout << std::setw(8) << "x" << std::setw(14) << "Aprox."
              << std::setw(14) << "Exato" << std::setw(14) << "Erro (%)\n";
    std::cout << std::string(50, '-') << "\n";
    for (int i = 0; i < M; i++)
    {
        double x = (i + 1) * h;
        double ex = exactSolutionPVC1(x);
        double err = (std::fabs(ex) > 1e-15) ? std::fabs((approx[i] - ex) / ex) * 100.0 : 0.0;
        std::cout << std::setw(8) << x << std::setw(14) << approx[i]
                  << std::setw(14) << ex << std::setw(14) << err << "\n";
    }
}

static void printPVC1Comparison(const std::string &methodName,
                                const std::vector<double> &approxN4,
                                const std::vector<double> &approxN8)
{
    double xs[3] = {0.25, 0.50, 0.75};
    int idxN8[3] = {1, 3, 5};
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n=== " << methodName << " - Comparacao N=4 vs N=8 ===\n";
    std::cout << std::setw(8) << "x" << std::setw(14) << "N=4"
              << std::setw(14) << "N=8" << std::setw(14) << "Exato"
              << std::setw(12) << "ErrN4(%)" << std::setw(12) << "ErrN8(%)\n";
    std::cout << std::string(74, '-') << "\n";
    for (int k = 0; k < 3; k++)
    {
        double x = xs[k];
        double ex = exactSolutionPVC1(x);
        double a4 = approxN4[k];
        double a8 = approxN8[idxN8[k]];
        double e4 = std::fabs((a4 - ex) / ex) * 100.0;
        double e8 = std::fabs((a8 - ex) / ex) * 100.0;
        std::cout << std::setw(8) << x << std::setw(14) << a4
                  << std::setw(14) << a8 << std::setw(14) << ex
                  << std::setw(12) << e4 << std::setw(12) << e8 << "\n";
    }
}

static void runTaskPVC1(const std::string &methodName, bool useFEM)
{
    std::cout << "\n=== Tarefa: PVC1 com N=8 (" << methodName << ") ===\n";
    auto approxN4 = useFEM ? finiteElementsPVC1(4) : finiteDiffPVC1(4);
    auto approxN8 = useFEM ? finiteElementsPVC1(8) : finiteDiffPVC1(8);
    printPVC1Table(methodName + " N=8", 8, approxN8);
    printPVC1Comparison(methodName, approxN4, approxN8);
}

static void runTaskPVC2_FD()
{
    std::cout << "\n=== Tarefa 19.2: PVC2 com N=8 (Diferencas Finitas) ===\n";
    std::cout << "PVC2: d2u/dx2 + d2u/dy2 = 4, u=0 no contorno\n";

    int N = 8;
    auto u8 = finiteDiffPVC2(N);
    auto u4 = finiteDiffPVC2(4);

    int M8 = N - 1;
    double h8 = 1.0 / N;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nResultados N=8:\n";
    std::cout << std::setw(6) << "j\\i";
    for (int i = 1; i <= M8; i++)
        std::cout << std::setw(12) << i * h8;
    std::cout << "\n"
              << std::string(6 + M8 * 12, '-') << "\n";

    for (int j = M8 - 1; j >= 0; j--)
    {
        std::cout << std::setw(6) << (j + 1) * h8;
        for (int i = 0; i < M8; i++)
            std::cout << std::setw(12) << u8[j * M8 + i];
        std::cout << "\n";
    }

    std::cout << "\n=== Comparacao N=4 vs N=8 nos mesmos pontos ===\n";
    std::cout << std::setw(8) << "x" << std::setw(8) << "y"
              << std::setw(14) << "N=4" << std::setw(14) << "N=8"
              << std::setw(14) << "Dif. relativa(%)\n";
    std::cout << std::string(58, '-') << "\n";

    int M4 = 3;
    double h4 = 0.25;
    int idx8[3] = {1, 3, 5};

    for (int jj = 0; jj < 3; jj++)
    {
        for (int ii = 0; ii < 3; ii++)
        {
            double x = (ii + 1) * h4, y = (jj + 1) * h4;
            double v4 = u4[jj * M4 + ii];
            double v8 = u8[idx8[jj] * M8 + idx8[ii]];
            double dr = (std::fabs(v4) > 1e-15) ? std::fabs((v8 - v4) / v4) * 100.0 : 0.0;
            std::cout << std::setw(8) << x << std::setw(8) << y
                      << std::setw(14) << v4 << std::setw(14) << v8
                      << std::setw(14) << dr << "\n";
        }
    }
}

// ---- sub-menus das aulas ----

static void menuFiniteDiff()
{
    while (true)
    {
        std::cout << "   Diferencas Finitas                     \n";
        std::cout << "  1 - PVC1: resolver com N escolhido      \n";
        std::cout << "  2 - Tarefa 19.1: PVC1 N=4 vs N=8        \n";
        std::cout << "  3 - Tarefa 19.2: PVC2 N=4 vs N=8        \n";
        std::cout << "  0 - Voltar                              \n";

        int op = readInt("Opcao: ", 0, 3);
        if (op == 0)
            return;
        if (op == 1)
        {
            int N = readInt("Numero de particoes N: ", 2, 1000);
            auto approx = finiteDiffPVC1(N);
            printPVC1Table("Diferencas Finitas", N, approx);
            pause();
        }
        else if (op == 2)
        {
            runTaskPVC1("Diferencas Finitas", false);
            pause();
        }
        else if (op == 3)
        {
            runTaskPVC2_FD();
            pause();
        }
    }
}

static void menuFiniteElements()
{
    while (true)
    {
        std::cout << "   Elementos Finitos                      \n";
        std::cout << "  1 - PVC1: resolver com N escolhido      \n";
        std::cout << "  2 - Tarefa 20.1: PVC1 N=4 vs N=8        \n";
        std::cout << "  0 - Voltar                              \n";

        int op = readInt("Opcao: ", 0, 2);
        if (op == 0)
            return;
        if (op == 1)
        {
            int N = readInt("Numero de elementos N: ", 2, 1000);
            auto approx = finiteElementsPVC1(N);
            printPVC1Table("Elementos Finitos", N, approx);
            pause();
        }
        else if (op == 2)
        {
            runTaskPVC1("Elementos Finitos", true);
            pause();
        }
    }
}

void menuPVC()
{
    while (true)
    {
        std::cout << "   PVC - Problemas de Valor de Contorno   \n";
        std::cout << "  1 - Diferencas Finitas                  \n";
        std::cout << "  2 - Elementos Finitos                   \n";
        std::cout << "  0 - Voltar                              \n";

        int op = readInt("Opcao: ", 0, 2);
        if (op == 0)
            return;
        if (op == 1)
            menuFiniteDiff();
        if (op == 2)
            menuFiniteElements();
    }
}