#include "../include/function.h"
#include "../include/diffLimits.h"

// ================ MÉTODOS DE DERIVAÇÃO INDIVIDUAIS ================
// -------- FORWARD, BACKWARD E CENTRAL --------

// FORWARD
float firstForwardDerivation(float x, float h){
    return ( F(x+h) - F(x) ) / h;
}

float secondForwardDerivation(float x, float h){
    return ( F(x + 2*h) - 2*F(x + h) + F(x) ) / (pow(h, 2));
}
    
float thirdFowardDerivation(float x, float h){
    return ( (1 / pow(h, 3)) * ( F(x+3*h) - 3*F(x+2*h) + 3*F(x+h) - F(x) ) );
}

// BACKWARD
float firstBackwardDerivation(float x, float h){
    return ( F(x) - F(x-h) ) / h;
}
    
float secondBackwardDerivation(float x, float h){
    return ( F(x) - F(x-h) ) / h;
}

// CENTRAL
float firstCentralDerivation(float x, float h){
    return ( F(x+h) - F(x-h) ) / (2*h);
}

float secondCentralDerivation(float x, float h){
    return ( F(x+h) - 2*F(x) + F(x-h) ) / (pow(h, 2));
}
// ================ MÉTODOS DE DERIVAÇÃO RECURSIVOS ================

// -------- FORWARD, BACKWARD E CENTRAL --------
float forwardRec(float x, float h, int order){
    if(order == 0) return F(x);

    float left = forwardRec(x+h, h, order - 1);
    float right = forwardRec(x, h, order - 1);

    return (left - right) / h;
}

float backwardRec(float x, float h, int order){
    if(order == 0) return F(x);

    float left = backwardRec(x, h, order - 1);
    float right = backwardRec(x-h, h, order - 1);

    return (left - right) / h;

}

float centralRec(float x, float h, int order){
    if (order == 0) return F(x);

    float left  = centralRec(x+h, h, order - 1);
    float right = centralRec(x-h, h, order - 1);

    return (left - right) / (2*h);
}