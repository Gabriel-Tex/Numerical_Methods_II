#include <iostream>
#include <string>
#include <limits>

#include "include/numericalDiff/menuNumericalDiff.h"
#include "include/numericalInt/menuNumericalInt.h"
#include "include/eigenvalues/menuEigenvalues.h"
#include "include/pvi/menuPVI.h"
#include "include/pvc/menuPVC.h"

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

int main()
{
    while (true)
    {
        std::cout << "  1 - Derivacao Numerica                  \n";
        std::cout << "  2 - Integracao Numerica                 \n";
        std::cout << "  3 - Autovalores e Autovetores           \n";
        std::cout << "  4 - PVI: Problemas de Valor Inicial     \n";
        std::cout << "  5 - PVC: Problemas de Valor de Contorno \n";
        std::cout << "  0 - Sair                                \n";

        int op = readInt("Opcao: ", 0, 5);
        switch (op)
        {
        case 1:
            menuNumericalDiff();
            break;
        case 2:
            menuNumericalInt();
            break;
        case 3:
            menuEigenvalues();
            break;
        case 4:
            menuPVI();
            break;
        case 5:
            menuPVC();
            break;
        case 0:
            std::cout << "Encerrando.\n";
            return 0;
        }
    }
}