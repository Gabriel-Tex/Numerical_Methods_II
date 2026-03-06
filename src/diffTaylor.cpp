#include "../include/function.h"
#include "../include/diffTaylor.h"
#include "../include/diffLimits.h"

// Fórmula de Taylor:
// F(x + h) = F(x) + F'(x)*h + 1/2! * F''(x)*h² + 1/3! * F'''(x)*h³ + ...

derivate firstTaylorDerivation_LinearErrorAndForward(double x, double h){

    derivate firstDerivate;
    firstDerivate.derivateValue = (F(x+h) - F(x)) / h;
    firstDerivate.error = -(1/2)*secondForwardDerivation(x, h)*h;

    return firstDerivate;
}

derivate firstTaylorDerivation_QuadraticErrorAndForward(double x, double h){
    // alfa * 1/2! * F''(x)*h² + beta * 1/2! * F''(x) * 2*h²
    // ( 1/2! * F''(x)*h² * [alfa + 4*beta]) => alfa + 4*beta = 0 => alfa = -4*beta
    // alfa = 4, beta = -1
    derivate firstDerivate;
    firstDerivate.derivateValue = ( 1/(2*h) ) * ( -F(x + 2*h) + 4*F(x+h)-3*F(x) );
    firstDerivate.error = (1/3) * thirdForwardDerivation(x, h) * pow(h, 2); // O(h²)

    return firstDerivate;
}

derivate firstTaylorDerivation_CubicErrorAndForward(double x, double h){
    // alfa = -1, beta = 1/2, gama = -1/9
    derivate firstDerivate;
    firstDerivate.derivateValue = ( (1/3) * F(x+3*h) - (3/2)*F(x+2*h) + 3*F(x + h) - 2*F(x) ) / h;

    firstDerivate.error = -(1/4) * fourthForwardDerivation(x, h) * pow(h, 3); // O(h³)

    return firstDerivate;
}

