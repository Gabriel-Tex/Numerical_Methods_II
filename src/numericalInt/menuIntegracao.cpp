#include "../../include/numericalInt/menuIntegracao.h"
#include "../../include/function.h"
#include "../../include/numericalInt/newtonCotes.h"
#include "../../include/numericalInt/gaussLegendre.h"
#include "../../include/numericalInt/gaussSpecial.h"
#include "../../include/numericalInt/singularIntegral.h"
#include "../../include/numericalInt/multiDimIntegral.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>

// ---- Utilitários locais ----

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
        std::cout << "  Opcao invalida. Tente novamente.\n";
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
        std::cout << "  Valor invalido. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// ================================ Newton-Cotes com particionamento ================================
void menuNewtonCotes()
{
    limpaTela();
    std::cout << "=== AULAS 7 e 8 - Newton-Cotes ===\n";
    std::cout << "Funcao: f(x) = (sin(2x) + 4x^2 + 3x)^2\n";
    std::cout << "Valor exato: I = 17.8764703\n\n";

    double a = lerDouble("Limite inferior a: ");
    double b = lerDouble("Limite superior b: ");

    std::cout << "\nAbordagem:\n  1 - Fechada\n  2 - Aberta\n";
    bool fechada = (lerInt("Escolha: ", 1, 2) == 1);

    int grau = lerInt("Grau do polinomio de substituicao (1-4): ", 1, 4);
    double tol = lerDouble("Tolerancia (ex: 1e-6): ");

    int iters = 0;
    double result = newtonCotesParticionado(F_NC, a, b, grau, fechada, tol, iters);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResultado: " << result << "\n";
    std::cout << "Particoes usadas (N): " << iters << "\n";
    pausar();
}

// ================================ Gauss-Legendre com particionamento ================================

void menuGaussLegendre()
{
    limpaTela();
    std::cout << "=== AULA 11 - Gauss-Legendre ===\n";
    std::cout << "Funcao: f(x) = (sin(2x) + 4x^2 + 3x)^2\n";
    std::cout << "Valor exato: I = 17.8764703\n\n";

    double a = lerDouble("Limite inferior a: ");
    double b = lerDouble("Limite superior b: ");

    int n = lerInt("Numero de pontos n (2, 3 ou 4): ", 2, 4);
    double tol = lerDouble("Tolerancia (ex: 1e-6): ");

    int iters = 0;
    double result = gaussLegendreParticionado(F_NC, a, b, n, tol, iters);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResultado: " << result << "\n";
    std::cout << "Particoes usadas (N): " << iters << "\n";
    pausar();
}

// ================================ Quadraturas Especiais de Gauss ================================
void menuGaussEspecial()
{
    limpaTela();
    std::cout << "=== AULA 12 - Quadraturas Especiais de Gauss ===\n\n";
    std::cout << "Quadratura:\n";
    std::cout << "  1 - Gauss-Hermite   (integral de -inf a +inf com peso e^{-s^2})\n";
    std::cout << "  2 - Gauss-Laguerre  (integral de 0 a +inf com peso e^{-s})\n";
    std::cout << "  3 - Gauss-Chebyshev (integral de -1 a +1 com peso 1/sqrt(1-s^2))\n";
    int quad = lerInt("Escolha: ", 1, 3);

    int n = lerInt("Numero de pontos n (2, 3 ou 4): ", 2, 4);

    std::cout << "\nUsando F(x) = x^3 (funcao padrao do projeto)\n";
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

// ================================ Integral Singular ================================
void menuSingular()
{
    limpaTela();
    std::cout << "=== AULA 14 - Integracao com Singularidade ===\n\n";
    std::cout << "Funcoes disponiveis:\n";
    std::cout << "  1 - f(x) = x^{-1/3}        (integral [0,1] = 3/2 = 1.5)\n";
    std::cout << "  2 - f(x) = x^{-1/2} e^{-x} (integral [0,+inf] = sqrt(pi) ~ 1.7724539)\n";
    int fn = lerInt("Escolha a funcao: ", 1, 2);

    std::cout << "\nMudanca de variavel:\n  1 - Exponencial simples\n  2 - Exponencial dupla\n";
    std::string tipoStr = (lerInt("Escolha: ", 1, 2) == 1) ? "simples" : "dupla";

    int nGL = lerInt("Pontos Gauss-Legendre (2, 3 ou 4): ", 2, 4);
    double cMin = lerDouble("Valor inicial de c (ex: 3.0): ");
    double tol = lerDouble("Tolerancia (ex: 1e-6): ");

    Func1D f;
    double a, b, exato;

    if (fn == 1)
    {
        f = [](double x)
        { return F_sing1(x); };
        a = 0.0;
        b = 1.0;
        exato = 1.5;
    }
    else
    {
        f = [](double x)
        { return F_sing2(x); };
        a = 0.0;
        b = 1.0;
        exato = std::sqrt(M_PI);
    }

    double cFinal = 0.0;
    double result = singularSolucao2(f, a, b, tipoStr, cMin, tol, nGL, cFinal);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\nResultado:     " << result << "\n";
    std::cout << "c final usado: " << cFinal << "\n";
    std::cout << "Valor exato:   " << exato << "\n";
    std::cout << "Erro absoluto: " << std::fabs(result - exato) << "\n";
    pausar();
}

// ================================Área de Superfície 3D ================================
void menuAreaSuperficie()
{
    limpaTela();
    std::cout << "=== AULA 15 - Area de Superficie 3D ===\n";
    std::cout << "Superficie: z = f(x,y) = 0.2*(x^2 - y^2)\n\n";

    std::cout << "Dominio:\n";
    std::cout << "  1 - Retangular: x em [-50,50], y em [-50,50]  (ref: 153467.00)\n";
    std::cout << "  2 - Eliptico:   x^2/1600 + y^2/1600 <= 1     (a=b=40)\n";
    int dom = lerInt("Escolha: ", 1, 2);

    int n = lerInt("Pontos GL em cada direcao (2, 3 ou 4): ", 2, 4);

    Func2D dfx = [](double x, double y)
    { return dF_surf_dx(x, y); };
    Func2D dfy = [](double x, double y)
    { return dF_surf_dy(x, y); };

    double result = 0.0;
    if (dom == 1)
    {
        result = areaSuperficieRetangular(dfx, dfy, -50.0, 50.0, -50.0, 50.0, n);
        std::cout << "\nResultado (dominio retangular): ";
    }
    else
    {
        result = areaSuperficieEliptica(dfx, dfy, 40.0, 40.0, n);
        std::cout << "\nResultado (dominio eliptico a=b=40): ";
    }

    std::cout << std::fixed << std::setprecision(4) << result << "\n";
    pausar();
}

// ================================ Volume abaixo de superfície ================================
void menuVolume()
{
    limpaTela();
    std::cout << "=== AULA 16 - Volume abaixo da Superficie ===\n";
    std::cout << "Superficie: z = f(x,y) = 0.2*(x^2 - y^2)\n\n";

    std::cout << "Dominio:\n";
    std::cout << "  1 - Retangular: x em [-40,40], y em [-20,20]  (ref: 256000 m^3)\n";
    std::cout << "  2 - Eliptico:   x^2/1600 + y^2/400 <= 1      (a=40, b=20)\n";
    int dom = lerInt("Escolha: ", 1, 2);

    int n = lerInt("Pontos GL em cada direcao (2, 3 ou 4): ", 2, 4);

    Func2D fsurf = [](double x, double y)
    { return F_surf(x, y); };

    double result = 0.0;
    if (dom == 1)
    {
        result = volumeRetangular(fsurf, -40.0, 40.0, -20.0, 20.0, n);
        std::cout << "\nResultado (dominio retangular): ";
    }
    else
    {
        result = volumeEliptico(fsurf, 40.0, 20.0, n);
        std::cout << "\nResultado (dominio eliptico a=40, b=20): ";
    }

    std::cout << std::fixed << std::setprecision(4) << result << " m^3\n";
    pausar();
}