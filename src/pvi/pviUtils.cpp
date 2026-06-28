#include "../../include/pvi/pviUtils.h"
#include <cmath>
#include <stdexcept>

// Soma elemento a elemento de dois estados
Estado somaEstado(const Estado &a, const Estado &b)
{
    int n = (int)a.size();
    Estado r(n);
    for (int i = 0; i < n; i++)
        r[i] = a[i] + b[i];
    return r;
}

// Multiplica todos os componentes do estado por um escalar
Estado escalaEstado(double alpha, const Estado &a)
{
    int n = (int)a.size();
    Estado r(n);
    for (int i = 0; i < n; i++)
        r[i] = alpha * a[i];
    return r;
}

// Erro relativo pela norma infinito: max_i |novo[i] - antigo[i]| / |novo[i]|
double erroRelativo(const Estado &novo, const Estado &antigo)
{
    double err = 0.0;
    int n = (int)novo.size();
    for (int i = 0; i < n; i++)
    {
        double den = std::fabs(novo[i]);
        if (den < 1e-15)
            den = 1e-15;
        double e = std::fabs(novo[i] - antigo[i]) / den;
        if (e > err)
            err = e;
    }
    return err;
}