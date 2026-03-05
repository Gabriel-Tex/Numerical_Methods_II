#ifndef DIFFLIMITS_H
#define DIFFLIMITS_H

// ================ MÉTODOS DE DERIVAÇÃO INDIVIDUAIS ================
// -------- FORWARD, BACKWARD E CENTRAL --------

// FORWARD
float firstForwardDerivation(float x, float h);

float secondForwardDerivation(float x, float h);
    
float thirdFowardDerivation(float x, float h);

// BACKWARD
float firstBackwardDerivation(float x, float h);
    
float secondBackwardDerivation(float x, float h);

// CENTRAL
float firstCentralDerivation(float x, float h);

float secondCentralDerivation(float x, float h);
// ================ MÉTODOS DE DERIVAÇÃO RECURSIVOS ================

// -------- FORWARD, BACKWARD E CENTRAL --------
float forwardRec(float x, float h, int order);

float backwardRec(float x, float h, int order);

float centralRec(float x, float h, int order);

#endif