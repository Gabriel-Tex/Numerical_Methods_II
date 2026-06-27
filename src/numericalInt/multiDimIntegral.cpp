#include "../../include/numericalInt/multiDimIntegral.h"
#include "../../include/numericalInt/gaussLegendre.h"
#include <cmath>

// ÁREA DE SUPERFÍCIE 3D - DOMÍNIO RETANGULAR
// A = ∫∫ √((∂f/∂x)² + (∂f/∂y)² + 1) dA
// Resolvido por GL 2D com n pontos em cada direção.
double areaSuperficieRetangular(Func2D dfx, Func2D dfy,
                                double xa, double xb,
                                double ya, double yb,
                                int n)
{
    // Integrando: g(x,y) = √((∂f/∂x)² + (∂f/∂y)² + 1)
    Func2D g = [&](double x, double y) -> double
    {
        double fx = dfx(x, y);
        double fy = dfy(x, y);
        return sqrt(fx * fx + fy * fy + 1.0);
    };
    return gaussLegendre2D(g, xa, xb, ya, yb, n);
}

// ÁREA DE SUPERFÍCIE 3D - DOMÍNIO ELÍPTICO
// Domínio: x²/a² + y²/b² ≤ 1
//
// Mudança de variável 1 (polar elíptica, como na elipse da Aula 15):
//   x(α,β) = a*α*cos(β),  y(α,β) = b*α*sin(β)
//   α ∈ [0,1],  β ∈ [0,2π]
//   |J₁| = a*b*α
//
// Mudança de variável 2 (GL → [-1,1]):
//   α(u) = (1+u)/2 → du = 2dα,  a_GL = [−1,1] → α = [0,1]
//   β(v) = π*(1+v) → dv = dβ/π, b_GL = [−1,1] → β = [0,2π]
//   |J₂| = (1/2) * π = π/2
//
// Integrando final (sobre [-1,1]²):
//   √(... + 1) * a*b*α * π/2
double areaSuperficieEliptica(Func2D dfx, Func2D dfy,
                              double a, double b, int n)
{
    Func2D integrand = [&](double u, double v) -> double
    {
        // Mudança de variável 2: u∈[-1,1]→α∈[0,1],  v∈[-1,1]→β∈[0,2π]
        double alpha = (1.0 + u) / 2.0;
        double beta = M_PI * (1.0 + v);

        // Coordenadas cartesianas via mudança 1
        double x = a * alpha * cos(beta);
        double y = b * alpha * sin(beta);

        // Jacobiano total: |J₁| * |J₂| = a*b*α * π/2
        double jac = a * b * alpha * (M_PI / 2.0);

        double fx = dfx(x, y);
        double fy = dfy(x, y);
        return sqrt(fx * fx + fy * fy + 1.0) * jac;
    };

    return gaussLegendre2D(integrand, -1.0, 1.0, -1.0, 1.0, n);
}

// VOLUME - DOMÍNIO RETANGULAR
// V = ∫∫ f(x,y) dA
double volumeRetangular(Func2D f,
                        double xa, double xb,
                        double ya, double yb,
                        int n)
{
    return gaussLegendre2D(f, xa, xb, ya, yb, n);
}

// VOLUME - DOMÍNIO ELÍPTICO
// Domínio: x²/a² + y²/b² ≤ 1
//
// Mudança 1: x=a*α*cos(β), y=b*α*sin(β), |J₁|=a*b*α
// Mudança 2: α=(1+u)/2, β=π*(1+v), |J₂|=π/2
double volumeEliptico(Func2D f, double a, double b, int n)
{
    Func2D integrand = [&](double u, double v) -> double
    {
        double alpha = (1.0 + u) / 2.0;
        double beta = M_PI * (1.0 + v);

        double x = a * alpha * cos(beta);
        double y = b * alpha * sin(beta);

        double jac = a * b * alpha * (M_PI / 2.0);

        return f(x, y) * jac;
    };

    return gaussLegendre2D(integrand, -1.0, 1.0, -1.0, 1.0, n);
}