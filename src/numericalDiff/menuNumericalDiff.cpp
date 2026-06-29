#include "../../include/numericalDiff/menuNumericalDiff.h"
#include "../../include/function.h"
#include "../../include/numericalDiff/finiteDiff.h"
#include "../../include/numericalDiff/finiteDiffNewton.h"
#include "../../include/numericalDiff/finiteDiffTaylor.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

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

void menuNumericalDiff()
{
    std::cout << "   Derivacao Numerica                     \n";
    std::cout << "  1 - Diferencas Finitas                  \n";
    std::cout << "  2 - Formula de Newton (2a deriv, 4a ord)\n";
    std::cout << "  3 - Expansao de Taylor (1a derivada)    \n";
    int metodo = readInt("Escolha: ", 1, 3);

    double x = readDouble("Ponto x: ");
    double h = readDouble("Passo h: ");

    std::cout << std::fixed << std::setprecision(10);

    if (metodo == 1)
    {
        std::cout << "\nTipo:\n  1 - Forward\n  2 - Backward\n  3 - Central\n";
        int tipo = readInt("Escolha: ", 1, 3);
        std::cout << "\nAbordagem:\n  1 - Formulas individuais\n  2 - Formula recursiva\n";
        int abord = readInt("Abordagem: ", 1, 2);
        int ordem = readInt("Ordem (1 a 4): ", 1, 4);

        double resultado = 0.0;
        if (abord == 1)
        {
            if (tipo == 1)
            {
                switch (ordem)
                {
                case 1:
                    resultado = firstForwardDerivation(x, h);
                    break;
                case 2:
                    resultado = secondForwardDerivation(x, h);
                    break;
                case 3:
                    resultado = thirdForwardDerivation(x, h);
                    break;
                case 4:
                    resultado = fourthForwardDerivation(x, h);
                    break;
                }
                std::cout << "\n[Forward, ordem " << ordem << ", formula individual]\n";
            }
            else if (tipo == 2)
            {
                switch (ordem)
                {
                case 1:
                    resultado = firstBackwardDerivation(x, h);
                    break;
                case 2:
                    resultado = secondBackwardDerivation(x, h);
                    break;
                case 3:
                    resultado = thirdBackwardDerivation(x, h);
                    break;
                case 4:
                    resultado = fourthBackwardDerivation(x, h);
                    break;
                }
                std::cout << "\n[Backward, ordem " << ordem << ", formula individual]\n";
            }
            else
            {
                switch (ordem)
                {
                case 1:
                    resultado = firstCentralDerivation(x, h);
                    break;
                case 2:
                    resultado = secondCentralDerivation(x, h);
                    break;
                case 3:
                    resultado = thirdCentralDerivation(x, h);
                    break;
                case 4:
                    resultado = fourthCentralDerivation(x, h);
                    break;
                }
                std::cout << "\n[Central, ordem " << ordem << ", formula individual]\n";
            }
        }
        else
        {
            if (tipo == 1)
            {
                resultado = forwardRec(x, h, ordem);
                std::cout << "\n[Forward, ordem " << ordem << ", recursiva]\n";
            }
            else if (tipo == 2)
            {
                resultado = backwardRec(x, h, ordem);
                std::cout << "\n[Backward, ordem " << ordem << ", recursiva]\n";
            }
            else
            {
                resultado = centralRec(x, h, ordem);
                std::cout << "\n[Central, ordem " << ordem << ", recursiva]\n";
            }
        }

        std::cout << "F'(x) aprox = " << resultado << "\n";
    }
    else if (metodo == 2)
    {
        double resultado = secondNewtonDerivation_fourthOrder(x, h);
        std::cout << "\n[Newton - 2a derivada, 4a ordem, forward]\n";
        std::cout << "F''(x) aprox = " << resultado << "\n";
    }
    else
    {
        std::cout << "\nOrdem do erro:\n";
        std::cout << "  1 - O(h)   - erro linear     (2 pontos, forward)\n";
        std::cout << "  2 - O(h^2) - erro quadratico (3 pontos, forward)\n";
        std::cout << "  3 - O(h^3) - erro cubico     (4 pontos, forward)\n";
        int ordem = readInt("Escolha: ", 1, 3);

        derivate d;
        std::string nomeOrdem;
        switch (ordem)
        {
        case 1:
            d = firstTaylorDerivation_LinearErrorAndForward(x, h);
            nomeOrdem = "O(h)";
            break;
        case 2:
            d = firstTaylorDerivation_QuadraticErrorAndForward(x, h);
            nomeOrdem = "O(h^2)";
            break;
        case 3:
            d = firstTaylorDerivation_CubicErrorAndForward(x, h);
            nomeOrdem = "O(h^3)";
            break;
        }

        std::cout << "\n[Taylor - 1a derivada, erro " << nomeOrdem << ", forward]\n";
        std::cout << "F'(x) aprox = " << d.derivateValue << "\n";
        std::cout << "Erro estimado = " << d.error << "\n";
    }

    pause();
}