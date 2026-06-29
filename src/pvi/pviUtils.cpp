#include "../../include/pvi/pviUtils.h"
#include <cmath>

State addStates(const State &a, const State &b)
{
    int n = (int)a.size();
    State r(n);
    for (int i = 0; i < n; i++)
        r[i] = a[i] + b[i];
    return r;
}

State scaleState(double alpha, const State &a)
{
    int n = (int)a.size();
    State r(n);
    for (int i = 0; i < n; i++)
        r[i] = alpha * a[i];
    return r;
}

double relativeError(const State &newer, const State &older)
{
    double err = 0.0;
    int n = (int)newer.size();
    for (int i = 0; i < n; i++)
    {
        double den = std::fabs(newer[i]);
        if (den < 1e-15)
            den = 1e-15;
        double e = std::fabs(newer[i] - older[i]) / den;
        if (e > err)
            err = e;
    }
    return err;
}