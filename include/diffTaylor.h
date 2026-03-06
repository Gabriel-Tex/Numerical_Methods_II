#ifndef DIFFTAYLOR_H
#define DIFFTAYLOR_H

#include <cmath>

struct derivate
{
    double derivateValue;
    double error;
};

derivate firstTaylorDerivation_LinearErrorAndForward(double x, double h);

derivate firstTaylorDerivation_QuadraticErrorAndForward(double x, double h);

derivate firstTaylorDerivation_CubicErrorAndForward(double x, double h);

#endif