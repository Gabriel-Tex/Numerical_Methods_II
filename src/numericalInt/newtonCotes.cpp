#include "../../include/numericalInt/newtonCotes.h"
#include <cmath>
#include <stdexcept>

// =========================== FECHADAS ===========================

// Grau 1 - Regra do Trapézio
// I ≈ (Δx/2) * [f(xi) + f(xf)]
double trapezioFechado(Func1D f, double xi, double xf)
{
    double dx = xf - xi;
    return (dx / 2.0) * (f(xi) + f(xf));
}

// Grau 2 - Regra de Simpson 1/3
// I ≈ (h/3) * [f(xi) + 4*f(xi+h) + f(xf)],  h = Δx/2
double simpson13(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 2.0;
    return (h / 3.0) * (f(xi) + 4.0 * f(xi + h) + f(xf));
}

// Grau 3 - Regra de Simpson 3/8
// I ≈ (3h/8) * [f(xi) + 3f(xi+h) + 3f(xi+2h) + f(xf)],  h = Δx/3
double simpson38(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 3.0;
    return (3.0 * h / 8.0) * (f(xi) + 3.0 * f(xi + h) + 3.0 * f(xi + 2.0 * h) + f(xf));
}

// Grau 4 - Regra de Boole (fechada)
// I ≈ (2h/45) * [7f(x0) + 32f(x1) + 12f(x2) + 32f(x3) + 7f(x4)],  h = Δx/4
double fechadaGrau4(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 4.0;
    return (2.0 * h / 45.0) * (7.0 * f(xi) + 32.0 * f(xi + h) + 12.0 * f(xi + 2.0 * h) + 32.0 * f(xi + 3.0 * h) + 7.0 * f(xf));
}

// =========================== ABERTAS ===========================

// Grau 1 - Trapézio Aberto
// I ≈ (Δx/2) * [f(xi+h) + f(xi+2h)],  h = Δx/3
double trapezioAberto(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 3.0;
    return ((xf - xi) / 2.0) * (f(xi + h) + f(xi + 2.0 * h));
}

// Grau 2 - Regra de Milne
// I ≈ (4h/3) * [2f(x1) - f(x2) + 2f(x3)],  h = Δx/4
double milne(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 4.0;
    return (4.0 * h / 3.0) * (2.0 * f(xi + h) - f(xi + 2.0 * h) + 2.0 * f(xi + 3.0 * h));
}

// Grau 3 - Aberta
// I ≈ (5h/24) * [11f(x1) + f(x2) + f(x3) + 11f(x4)],  h = Δx/5
double abertaGrau3(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 5.0;
    return (5.0 * h / 24.0) * (11.0 * f(xi + h) + 1.0 * f(xi + 2.0 * h) + 1.0 * f(xi + 3.0 * h) + 11.0 * f(xi + 4.0 * h));
}

// Grau 4 - Aberta 
// Pontos internos: x1..x5,  h = Δx/6
// I ≈ (6h/20) * [11f(x1) - 14f(x2) + 26f(x3) - 14f(x4) + 11f(x5)]
double abertaGrau4(Func1D f, double xi, double xf)
{
    double h = (xf - xi) / 6.0;
    return (6.0 * h / 20.0) * (11.0 * f(xi + h) - 14.0 * f(xi + 2.0 * h) + 26.0 * f(xi + 3.0 * h) - 14.0 * f(xi + 4.0 * h) + 11.0 * f(xi + 5.0 * h));
}

// =========================== PARTICIONAMENTO ADAPTATIVO ===========================
// Duplica N a cada iteração até |I_novo - I_ant| < tolerancia

double newtonCotesParticionado(Func1D f, double a, double b, int grau, bool fechada,
                                    double tolerancia, int &iteracoes)
{
    auto formula = [&](Func1D fn, double xi, double xf) -> double
    {
        if (fechada)
        {
            switch (grau)
            {
            case 1:
                return trapezioFechado(fn, xi, xf);
            case 2:
                return simpson13(fn, xi, xf);
            case 3:
                return simpson38(fn, xi, xf);
            case 4:
                return fechadaGrau4(fn, xi, xf);
            default:
                throw std::invalid_argument("Grau inválido");
            }
        }
        else
        {
            switch (grau)
            {
            case 1:
                return trapezioAberto(fn, xi, xf);
            case 2:
                return milne(fn, xi, xf);
            case 3:
                return abertaGrau3(fn, xi, xf);
            case 4:
                return abertaGrau4(fn, xi, xf);
            default:
                throw std::invalid_argument("Grau inválido");
            }
        }
    };

    int N = 1;
    double dx = (b - a) / N;
    double I_ant = formula(f, a, b);

    for (int iter = 1; iter <= 100000; iter++)
    {
        N *= 2;
        dx = (b - a) / N;
        double I_novo = 0.0;
        for (int i = 0; i < N; i++)
        {
            double xi = a + i * dx;
            double xf = xi + dx;
            I_novo += formula(f, xi, xf);
        }
        iteracoes = N;
        if (fabs(I_novo - I_ant) < tolerancia)
        {
            return I_novo;
        }
        I_ant = I_novo;
    }
    return I_ant;
}