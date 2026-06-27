#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <limits>

#include "include/function.h"
#include "include/numericalInt/newtonCotes.h"
#include "include/numericalInt/gaussLegendre.h"
#include "include/numericalInt/gaussSpecial.h"
#include "include/numericalInt/singularIntegral.h"
#include "include/numericalInt/multiDimIntegral.h"

static void limpaTela()
{
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

static void pausar()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int lerInt(const std::string &prompt, int min, int max)
{
    int v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v && v >= min && v <= max)
            return v;
        std::cout << "  Opção inválida. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static double lerDouble(const std::string &prompt)
{
    double v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v)
            return v;
        std::cout << "  Valor inválido. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// AULAS 7 e 8 - Newton-Cotes com particionamento

static void menuNewtonCotes()
{
    limpaTela();
    std::cout << "=== AULAS 7 e 8 - Newton-Cotes ===\n";
    std::cout << "Função: f(x) = (sin(2x) + 4x² + 3x)²\n";
    std::cout << "Valor exato: I = 17.8764703\n\n";

    double a = lerDouble("Limite inferior a: ");
    double b = lerDouble("Limite superior b: ");

    std::cout << "\nAbordagem:\n  1 - Fechada\n  2 - Aberta\n";
    int abord = lerInt("Escolha: ", 1, 2);
    bool fechada = (abord == 1);

    std::cout << "\nGrau do polinômio de substituição (1-4): ";
    int grau = lerInt("Grau: ", 1, 4);

    double tol = lerDouble("Tolerância (ex: 1e-6): ");

    int iters = 0;
    double result = newtonCotesParticionado(F_NC, a, b, grau, fechada, tol, iters);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResultado: " << result << "\n";
    std::cout << "Partições usadas (N): " << iters << "\n";
    pausar();
}

// AULA 11 - Gauss-Legendre com particionamento

static void menuGaussLegendre()
{
    limpaTela();
    std::cout << "=== AULA 11 - Gauss-Legendre ===\n";
    std::cout << "Função: f(x) = (sin(2x) + 4x² + 3x)²\n";
    std::cout << "Valor exato: I = 17.8764703\n\n";

    double a = lerDouble("Limite inferior a: ");
    double b = lerDouble("Limite superior b: ");

    int n = lerInt("Número de pontos n (2, 3 ou 4): ", 2, 4);
    double tol = lerDouble("Tolerância (ex: 1e-6): ");

    int iters = 0;
    double result = gaussLegendreParticionado(F_NC, a, b, n, tol, iters);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResultado: " << result << "\n";
    std::cout << "Partições usadas (N): " << iters << "\n";
    pausar();
}

// AULA 12 - Quadraturas especiais (Hermite, Laguerre, Chebyshev)

static void menuGaussEspecial()
{
    limpaTela();
    std::cout << "=== AULA 12 - Quadraturas Especiais de Gauss ===\n\n";
    std::cout << "Quadratura:\n";
    std::cout << "  1 - Gauss-Hermite   (∫ e^{-s²} f(s) ds, limites: -∞ a +∞)\n";
    std::cout << "  2 - Gauss-Laguerre  (∫ e^{-s}  f(s) ds, limites: 0 a +∞)\n";
    std::cout << "  3 - Gauss-Chebyshev (∫ 1/√(1-s²) f(s) ds, limites: -1 a +1)\n";
    int quad = lerInt("Escolha: ", 1, 3);

    int n = lerInt("Número de pontos n (2, 3 ou 4): ", 2, 4);

    std::cout << "\nUsando F(x) = x³ (função padrão do projeto)\n";
    Func1D f = [](double s)
    { return F(s); };

    double result = 0.0;
    std::string nome;
    switch (quad)
    {
    case 1:
        result = gaussHermite(f, n);
        nome = "Gauss-Hermite";
        break;
    case 2:
        result = gaussLaguerre(f, n);
        nome = "Gauss-Laguerre";
        break;
    case 3:
        result = gaussChebyshev(f, n);
        nome = "Gauss-Chebyshev";
        break;
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n"
              << nome << " (n=" << n << "): " << result << "\n";
    pausar();
}

// AULA 14 - Integral singular (Solução 2: exp simples ou dupla)

static void menuSingular()
{
    limpaTela();
    std::cout << "=== AULA 14 - Integração com Singularidade ===\n\n";
    std::cout << "Funções disponíveis:\n";
    std::cout << "  1 - f(x) = x^{-1/3}      (integral [0,1] = 3/2 = 1.5)\n";
    std::cout << "  2 - f(x) = x^{-1/2}e^{-x} (integral [0,+∞] = √π ≈ 1.7724539)\n";
    int fn = lerInt("Escolha a função: ", 1, 2);

    std::cout << "\nMudança de variável:\n  1 - Exponencial simples\n  2 - Exponencial dupla\n";
    int tipo = lerInt("Escolha: ", 1, 2);
    std::string tipoStr = (tipo == 1) ? "simples" : "dupla";

    int nGL = lerInt("Pontos Gauss-Legendre (2, 3 ou 4): ", 2, 4);
    double cMin = lerDouble("Valor inicial de c (ex: 3.0): ");
    double tol = lerDouble("Tolerância (ex: 1e-6): ");

    Func1D f;
    double a, b;
    double exato;

    if (fn == 1)
    {
        // Singularidade em x=0, integral de 0 a 1
        f = [](double x)
        { return F_sing1(x); };
        a = 0.0;
        b = 1.0;
        exato = 1.5;
    }
    else
    {
        // Singularidade em x=0, integral de 0 a +∞
        // Usamos exp dupla com a=0, b=1 (o mapeamento leva [0,1]→[0,+∞] de forma suave)
        // Para integral [0,+∞]: usamos b muito grande e exp dupla
        f = [](double x)
        { return F_sing2(x); };
        a = 0.0;
        b = 1.0; // serão expandidos pela mudança de variável
        exato = sqrt(M_PI);
    }

    double cFinal = 0.0;
    double result = singularSolucao2(f, a, b, tipoStr, cMin, tol, nGL, cFinal);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\nResultado: " << result << "\n";
    std::cout << "c final usado: " << cFinal << "\n";
    std::cout << "Valor exato:  " << exato << "\n";
    std::cout << "Erro absoluto: " << fabs(result - exato) << "\n";
    pausar();
}

// AULA 15 - Área de superfície 3D

static void menuAreaSuperficie()
{
    limpaTela();
    std::cout << "=== AULA 15 - Área de Superfície 3D ===\n";
    std::cout << "Superfície: z = f(x,y) = 0.2*(x² - y²)\n\n";

    std::cout << "Domínio:\n";
    std::cout << "  1 - Retangular: x∈[-50,50], y∈[-50,50]  (Problema 1, ref: 153467.00)\n";
    std::cout << "  2 - Elíptico:   x²/1600 + y²/1600 ≤ 1   (Tarefa Aula 15, a=b=40)\n";
    int dom = lerInt("Escolha: ", 1, 2);

    int n = lerInt("Pontos GL em cada direção (2, 3 ou 4): ", 2, 4);

    Func2D dfx = [](double x, double y)
    { return dF_surf_dx(x, y); };
    Func2D dfy = [](double x, double y)
    { return dF_surf_dy(x, y); };

    double result = 0.0;
    if (dom == 1)
    {
        result = areaSuperficieRetangular(dfx, dfy, -50.0, 50.0, -50.0, 50.0, n);
        std::cout << "\nResultado (domínio retangular): ";
    }
    else
    {
        double a = 40.0, b = 40.0;
        result = areaSuperficieEliptica(dfx, dfy, a, b, n);
        std::cout << "\nResultado (domínio elíptico a=b=40): ";
    }

    std::cout << std::fixed << std::setprecision(4) << result << "\n";
    pausar();
}

// AULA 16 - Volume abaixo de superfície

static void menuVolume()
{
    limpaTela();
    std::cout << "=== AULA 16 - Volume abaixo da Superfície ===\n";
    std::cout << "Superfície: z = f(x,y) = 0.2*(x² - y²)\n\n";

    std::cout << "Domínio:\n";
    std::cout << "  1 - Retangular: x∈[-40,40], y∈[-20,20]  (Problema 1, ref: 256000 m³)\n";
    std::cout << "  2 - Elíptico:   x²/1600 + y²/400 ≤ 1    (Tarefa Aula 16, a=40, b=20)\n";
    int dom = lerInt("Escolha: ", 1, 2);

    int n = lerInt("Pontos GL em cada direção (2, 3 ou 4): ", 2, 4);

    Func2D fsurf = [](double x, double y)
    { return F_surf(x, y); };

    double result = 0.0;
    if (dom == 1)
    {
        result = volumeRetangular(fsurf, -40.0, 40.0, -20.0, 20.0, n);
        std::cout << "\nResultado (domínio retangular): ";
    }
    else
    {
        double a = 40.0, b = 20.0;
        result = volumeEliptico(fsurf, a, b, n);
        std::cout << "\nResultado (domínio elíptico a=40, b=20): ";
    }

    std::cout << std::fixed << std::setprecision(4) << result << " m³\n";
    pausar();
}

// Menu principal

int main()
{
    while (true)
    {
        limpaTela();
        std::cout << "╔══════════════════════════════════════════════════╗\n";
        std::cout << "║         MÉTODOS NUMÉRICOS II - INTEGRAÇÃO        ║\n";
        std::cout << "╠══════════════════════════════════════════════════╣\n";
        std::cout << "║  1 - Aulas  7/8  Newton-Cotes (com partição)     ║\n";
        std::cout << "║  2 - Aula   11   Gauss-Legendre (com partição)   ║\n";
        std::cout << "║  3 - Aula   12   Quadraturas Especiais de Gauss  ║\n";
        std::cout << "║  4 - Aula   14   Integral Singular (exp. simples/║\n";
        std::cout << "║                  dupla - Solução 2)              ║\n";
        std::cout << "║  5 - Aula   15   Área de Superfície 3D           ║\n";
        std::cout << "║  6 - Aula   16   Volume                          ║\n";
        std::cout << "║  0 - Sair                                        ║\n";
        std::cout << "╚══════════════════════════════════════════════════╝\n";

        int op = lerInt("Opção: ", 0, 6);
        switch (op)
        {
        case 1:
            menuNewtonCotes();
            break;
        case 2:
            menuGaussLegendre();
            break;
        case 3:
            menuGaussEspecial();
            break;
        case 4:
            menuSingular();
            break;
        case 5:
            menuAreaSuperficie();
            break;
        case 6:
            menuVolume();
            break;
        case 0:
            std::cout << "Encerrando.\n";
            return 0;
        }
    }
}