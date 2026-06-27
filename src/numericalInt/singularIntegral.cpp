#include "../../include/numericalInt/singularIntegral.h"
#include "../../include/numericalInt/gaussLegendre.h"
#include <cmath>
#include <string>

// ========================== Mudança de variável - exponencial simples ==========================
double xSimples(double s, double a, double b)
{
    return (a + b) / 2.0 + (b - a) / 2.0 * tanh(s);
}

double dxSimples(double s, double a, double b)
{
    double cs = cosh(s);
    return (b - a) / 2.0 * (1.0 / (cs * cs)); // sech²(s)
}

// ========================== Mudança de variável - exponencial dupla ==========================
double xDupla(double s, double a, double b)
{
    return (a + b) / 2.0 + (b - a) / 2.0 * tanh(M_PI / 2.0 * sinh(s));
}

double dxDupla(double s, double a, double b)
{
    double cs = cosh(s);
    double arg = M_PI / 2.0 * sinh(s);
    double cd = cosh(arg);
    return (b - a) / 2.0 * (M_PI / 2.0 * cs / (cd * cd));
}

// Aproxima ∫_{-∞}^{+∞} f_bar(s) ds por ∫_{-c}^{+c} f_bar(s) ds
// f_bar(s) = f(x(s)) * dx/ds
// Usa Gauss-Legendre com nGL pontos e aumenta c até convergir.
double singularSolucao2(Func1D f, double a, double b,
                        const std::string &tipo,
                        double cMin, double tolerancia,
                        int nGL, double &cFinal)
{

    // Monta o integrando transformado f_bar(s)
    Func1D fBar;
    if (tipo == "simples")
    {
        fBar = [&](double s) -> double
        {
            double x = xSimples(s, a, b);
            return f(x) * dxSimples(s, a, b);
        };
    }
    else
    {
        fBar = [&](double s) -> double
        {
            // Proteção contra overflow em cosh(π/2 * sinh(s))
            double arg = M_PI / 2.0 * sinh(s);
            if (fabs(arg) > 700.0)
                return 0.0;
            double x = xDupla(s, a, b);
            double dx = dxDupla(s, a, b);
            if (!std::isfinite(x) || !std::isfinite(dx))
                return 0.0;
            double fx = f(x);
            if (!std::isfinite(fx))
                return 0.0;
            return fx * dx;
        };
    }

    double c = cMin;
    double I_ant = gaussLegendre(fBar, -c, +c, nGL);

    for (int iter = 0; iter < 1000; iter++)
    {
        c += 1.0;
        double I_novo = gaussLegendre(fBar, -c, +c, nGL);
        if (fabs(I_novo - I_ant) < tolerancia)
        {
            cFinal = c;
            return I_novo;
        }
        I_ant = I_novo;
    }
    cFinal = c;
    return I_ant;
}