#include "../../include/pvc/menuPVC.h"
#include "../../include/pvc/diferencasFinitas.h"
#include "../../include/pvc/elementosFinitos.h"
#include "../../include/pvc/pvcUtils.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// ---------- utilitarios locais ----------

static void pausarPVC()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int lerIntPVC(const std::string &prompt, int min, int max)
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

// ---------- exibicao de tabelas ----------

// Tabela comparativa para PVC1 (uma coluna de metodo)
static void exibeTabelaPVC1(
    const std::string &nomeMetodo,
    int N,
    const std::vector<double> &aprox)
{
    int M = N - 1;
    double h = 1.0 / N;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n--- " << nomeMetodo << " (N=" << N << ") ---\n";
    std::cout << std::setw(8)  << "x"
              << std::setw(14) << "Aprox."
              << std::setw(14) << "Exato"
              << std::setw(14) << "Erro (%)\n";
    std::cout << std::string(50, '-') << "\n";

    for (int i = 0; i < M; i++)
    {
        double x = (i + 1) * h;
        double ex = solucaoExataPVC1(x);
        double err = (std::fabs(ex) > 1e-15)
                     ? std::fabs((aprox[i] - ex) / ex) * 100.0
                     : 0.0;
        std::cout << std::setw(8)  << x
                  << std::setw(14) << aprox[i]
                  << std::setw(14) << ex
                  << std::setw(14) << err << "\n";
    }
}

// Exibe tabela ampliada conforme solicitado nas Tarefas 19 e 20
// Compara N=4 (da aula) e N=8 (tarefa)
static void exibeTabelaComparadaPVC1(
    const std::string &nomeMetodo,
    std::vector<double> aproxN4,
    std::vector<double> aproxN8)
{
    // Nos em comum: x = 0.25, 0.5, 0.75 (indices 0,1,2 para N=4)
    //              x = 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875 (N=8)
    // A tabela exibe apenas os nos do N=4 tambem presentes em N=8
    // (x = 0.25 -> indice 1 no N=8, x=0.5 -> indice 3, x=0.75 -> indice 5)

    double xs[3]      = {0.25, 0.50, 0.75};
    int    idxN8[3]   = {1, 3, 5};

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n=== " << nomeMetodo << " - Comparacao N=4 vs N=8 ===\n";
    std::cout << std::setw(8)  << "x"
              << std::setw(14) << "N=4"
              << std::setw(14) << "N=8"
              << std::setw(14) << "Exato"
              << std::setw(12) << "ErrN4(%)"
              << std::setw(12) << "ErrN8(%)\n";
    std::cout << std::string(74, '-') << "\n";

    for (int k = 0; k < 3; k++)
    {
        double x  = xs[k];
        double ex = solucaoExataPVC1(x);
        double a4 = aproxN4[k];
        double a8 = aproxN8[idxN8[k]];
        double e4 = std::fabs((a4 - ex) / ex) * 100.0;
        double e8 = std::fabs((a8 - ex) / ex) * 100.0;

        std::cout << std::setw(8)  << x
                  << std::setw(14) << a4
                  << std::setw(14) << a8
                  << std::setw(14) << ex
                  << std::setw(12) << e4
                  << std::setw(12) << e8 << "\n";
    }
}

// ---------- menus das tarefas ----------

// Tarefa 19.1 e 20.1: PVC1 com N=8
static void menuTarefa19_20_PVC1(const std::string &nomeMetodo, bool mef)
{
    std::cout << "\n=== Tarefa: PVC1 com N=8 (" << nomeMetodo << ") ===\n";

    // N=4 (resultado da aula, para comparacao)
    auto aproxN4 = mef ? elementosFinitasPVC1(4) : diferencasFinitasPVC1(4);
    // N=8 (tarefa)
    auto aproxN8 = mef ? elementosFinitasPVC1(8) : diferencasFinitasPVC1(8);

    exibeTabelaPVC1(nomeMetodo + " N=8", 8, aproxN8);
    exibeTabelaComparadaPVC1(nomeMetodo, aproxN4, aproxN8);
}

// Tarefa 19.2: PVC2 com N=8 (apenas Diferencas Finitas)
static void menuTarefa19_PVC2()
{
    std::cout << "\n=== Tarefa 19.2: PVC2 com N=8 (Diferencas Finitas) ===\n";
    std::cout << "PVC2: d2u/dx2 + d2u/dy2 = 4, u=0 no contorno\n";

    int N = 8;
    auto u8 = diferencasFinitasPVC2(N);

    // Referencia: N=4 (da aula)
    auto u4 = diferencasFinitasPVC2(4);

    // Exibe todos os nos interiores de N=8
    int M8 = N - 1;
    double h8 = 1.0 / N;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nResultados N=8:\n";
    std::cout << std::setw(6) << "j\\i";
    for (int i = 1; i <= M8; i++)
        std::cout << std::setw(12) << i * h8;
    std::cout << "\n" << std::string(6 + M8 * 12, '-') << "\n";

    for (int j = M8 - 1; j >= 0; j--)
    {
        std::cout << std::setw(6) << (j + 1) * h8;
        for (int i = 0; i < M8; i++)
            std::cout << std::setw(12) << u8[j * M8 + i];
        std::cout << "\n";
    }

    // Comparacao nos pontos do N=4 presentes em N=8
    // Para N=4: nos em x=0.25, 0.5, 0.75 e y=0.25, 0.5, 0.75 (3x3 = 9 nos)
    // Para N=8: mesmos nos tem indices i=1,3,5 e j=1,3,5
    std::cout << "\n=== Comparacao N=4 vs N=8 nos mesmos pontos ===\n";
    std::cout << std::setw(8)  << "x"
              << std::setw(8)  << "y"
              << std::setw(14) << "N=4"
              << std::setw(14) << "N=8"
              << std::setw(14) << "Dif. relativa(%)\n";
    std::cout << std::string(58, '-') << "\n";

    int M4 = 3;
    double h4 = 0.25;
    int idx8[3] = {1, 3, 5};

    for (int jj = 0; jj < 3; jj++)
    {
        for (int ii = 0; ii < 3; ii++)
        {
            double x = (ii + 1) * h4;
            double y = (jj + 1) * h4;
            double v4 = u4[jj * M4 + ii];
            double v8 = u8[idx8[jj] * M8 + idx8[ii]];
            double dr = (std::fabs(v4) > 1e-15)
                        ? std::fabs((v8 - v4) / v4) * 100.0 : 0.0;
            std::cout << std::setw(8) << x
                      << std::setw(8) << y
                      << std::setw(14) << v4
                      << std::setw(14) << v8
                      << std::setw(14) << dr << "\n";
        }
    }
}

// ---------- sub-menus das aulas ----------

// Aula 27: Diferencas Finitas
static void menuAula27()
{
    while (true)
    {
        std::cout << "\n╔══════════════════════════════════════════╗\n";
        std::cout << "║   Diferencas Finitas (Aula 27)           ║\n";
        std::cout << "╠══════════════════════════════════════════╣\n";
        std::cout << "║  1 - PVC1: resolver com N escolhido      ║\n";
        std::cout << "║  2 - Tarefa 19.1: PVC1 N=4 vs N=8       ║\n";
        std::cout << "║  3 - Tarefa 19.2: PVC2 N=4 vs N=8       ║\n";
        std::cout << "║  0 - Voltar                              ║\n";
        std::cout << "╚══════════════════════════════════════════╝\n";

        int op = lerIntPVC("Opcao: ", 0, 3);
        if (op == 0) return;

        if (op == 1)
        {
            int N = lerIntPVC("Numero de particoes N: ", 2, 1000);
            auto aprox = diferencasFinitasPVC1(N);
            exibeTabelaPVC1("Diferencas Finitas", N, aprox);
            pausarPVC();
        }
        else if (op == 2)
        {
            menuTarefa19_20_PVC1("Diferencas Finitas", false);
            pausarPVC();
        }
        else if (op == 3)
        {
            menuTarefa19_PVC2();
            pausarPVC();
        }
    }
}

// Aula 28: Elementos Finitos
static void menuAula28()
{
    while (true)
    {
        std::cout << "\n╔══════════════════════════════════════════╗\n";
        std::cout << "║   Elementos Finitos (Aula 28)            ║\n";
        std::cout << "╠══════════════════════════════════════════╣\n";
        std::cout << "║  1 - PVC1: resolver com N escolhido      ║\n";
        std::cout << "║  2 - Tarefa 20.1: PVC1 N=4 vs N=8       ║\n";
        std::cout << "║  0 - Voltar                              ║\n";
        std::cout << "╚══════════════════════════════════════════╝\n";

        int op = lerIntPVC("Opcao: ", 0, 2);
        if (op == 0) return;

        if (op == 1)
        {
            int N = lerIntPVC("Numero de elementos N: ", 2, 1000);
            auto aprox = elementosFinitasPVC1(N);
            exibeTabelaPVC1("Elementos Finitos", N, aprox);
            pausarPVC();
        }
        else if (op == 2)
        {
            menuTarefa19_20_PVC1("Elementos Finitos", true);
            pausarPVC();
        }
    }
}

// ---------- menu principal do PVC ----------

void menuPVC()
{
    while (true)
    {
        std::cout << "\n╔══════════════════════════════════════════╗\n";
        std::cout << "║   PVC - Problemas de Valor de Contorno   ║\n";
        std::cout << "╠══════════════════════════════════════════╣\n";
        std::cout << "║  1 - Diferencas Finitas (Aula 27)        ║\n";
        std::cout << "║  2 - Elementos Finitos  (Aula 28)        ║\n";
        std::cout << "║  0 - Voltar                              ║\n";
        std::cout << "╚══════════════════════════════════════════╝\n";

        int op = lerIntPVC("Opcao: ", 0, 2);
        if (op == 0) return;
        if (op == 1) menuAula27();
        if (op == 2) menuAula28();
    }
}