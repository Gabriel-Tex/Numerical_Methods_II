#include <iostream>
#include <string>
#include <limits>

#include "include/numericalDiff/menuDerivacao.h"
#include "include/numericalInt/menuIntegracao.h"
#include "include/eigenValues/menuAutovalores.h"

static void limpaTela()
{
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

static int lerInt(const std::string &prompt, int min, int max)
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
        limpaTela();
        std::cout << "══════════════════════════════════════\n";
        std::cout << "  1 - Derivacao Numerica              \n";
        std::cout << "  2 - Newton-Cotes (com particao)     \n";
        std::cout << "  3 - Gauss-Legendre (com particao)   \n";
        std::cout << "  4 - Quadraturas Especiais de Gauss  \n";
        std::cout << "  5 - Integral Singular               \n";
        std::cout << "  6 - Area de Superficie 3D           \n";
        std::cout << "  7 - Volume                          \n";
        std::cout << "──────────────────────────────────────\n";
        std::cout << "  8 - Metodos da Potencia             \n";
        std::cout << "  9 - Householder (Tarefa 13)         \n";
        std::cout << "  10 - Jacobi (Tarefa 14)             \n";
        std::cout << "  11 - QR (Tarefa 15)                 \n";
        std::cout << "  0 - Sair                            \n";
        std::cout << "══════════════════════════════════════\n";

        int op = lerInt("Opcao: ", 0, 10);
        switch (op)
        {
        case 1:  
            menuDerivacao();   
            break;
        case 2:
            menuNewtonCotes();
            break;
        case 3:
            menuGaussLegendre();
            break;
        case 4:
            menuGaussEspecial();
            break;
        case 5:
            menuSingular();
            break;
        case 6:
            menuAreaSuperficie();
            break;
        case 7:
            menuVolume();
            break;
        case 8:
            menuPotencia();
            break;
        case 9:
            menuHouseholder();
            break;
        case 10:
            menuJacobi();
            break;
        case 11:
            menuQR();
            break;
        case 0:
            std::cout << "Encerrando.\n";
            return 0;
        }
    }
}