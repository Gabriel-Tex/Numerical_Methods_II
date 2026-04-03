#include "../../include/function.h"

// corrigir: 
double secondNewtonDerivation_fourthOrderAndFoward(double x, double h){
    return (-(5/6)*F(x+5*h) + (61/12)*F(x+4*h) - 13*F(x+3*h) + (107/6)*F(x+2*h) - 
    (77/6)*F(x+h) + (45/12)*F(x)) / (h*h);
}