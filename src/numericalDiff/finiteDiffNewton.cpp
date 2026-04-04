#include "../../include/function.h"

double secondNewtonDerivation_fourthOrder(double x, double h){
    return ( 137*F(x + 6*h) -672*F(x + 5*h) +1635*F(x + 4*h) -2740*F(x + 3*h)
+3255*F(x + 2*h) - 2292*F(x + h) + 677*F(x) ) / (180 * h*h);
}

