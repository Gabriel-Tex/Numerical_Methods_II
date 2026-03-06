#include "../include/function.h"
#include "../include/diffLimits.h"

// ================ MÉTODOS DE DERIVAÇÃO INDIVIDUAIS ================
// -------- FORWARD, BACKWARD E CENTRAL --------

// FORWARD
double firstForwardDerivation(double x, double h){
    return ( F(x+h) - F(x) ) / h;
}

double secondForwardDerivation(double x, double h){
    return ( F(x + 2*h) - (2*F(x + h)) + F(x) ) / (pow(h, 2));
}
    
double thirdForwardDerivation(double x, double h){
    return ( (1 / pow(h, 3)) * ( F(x+3*h) - 3*F(x+2*h) + 3*F(x+h) - F(x) ) );
}

// BACKWARD
double firstBackwardDerivation(double x, double h){
    return ( F(x) - F(x-h) ) / h;
}
    
double secondBackwardDerivation(double x, double h){
    return ( F(x) - 2*F(x-h) + F(x-2*h) ) / pow(h, 2);
}

// CENTRAL
double firstCentralDerivation(double x, double h){
    return ( F(x+h) - F(x-h) ) / (2*h);
}

double secondCentralDerivation(double x, double h){
    return ( F(x+h) - 2*F(x) + F(x-h) ) / (pow(h, 2));
}
// ================ MÉTODOS DE DERIVAÇÃO RECURSIVOS ================

// -------- FORWARD, BACKWARD E CENTRAL --------
double forwardRec(double x, double h, int order){
    if(order == 0) return F(x);

    double left = forwardRec(x+h, h, order - 1);
    double right = forwardRec(x, h, order - 1);

    return (left - right) / h;
}

double backwardRec(double x, double h, int order){
    if(order == 0) return F(x);

    double left = backwardRec(x, h, order - 1);
    double right = backwardRec(x-h, h, order - 1);

    return (left - right) / h;

}

double centralRec(double x, double h, int order){
    if (order == 0) return F(x);
    
    double left  = centralRec(x+h, h, order - 1);
    double right = centralRec(x-h, h, order - 1);

    return (left - right) / (2*h);
}