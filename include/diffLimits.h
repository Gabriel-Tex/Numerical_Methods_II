#ifndef DIFFLIMITS_H
#define DIFFLIMITS_H

#include <cmath>

// ================ MÉTODOS DE DERIVAÇÃO INDIVIDUAIS ================
// -------- FORWARD, BACKWARD E CENTRAL --------

// FORWARD
double firstForwardDerivation(double x, double h);

double secondForwardDerivation(double x, double h);
    
double thirdForwardDerivation(double x, double h);

// BACKWARD
double firstBackwardDerivation(double x, double h);
    
double secondBackwardDerivation(double x, double h);

// CENTRAL
double firstCentralDerivation(double x, double h);

double secondCentralDerivation(double x, double h);
// ================ MÉTODOS DE DERIVAÇÃO RECURSIVOS ================

// -------- FORWARD, BACKWARD E CENTRAL --------
double forwardRec(double x, double h, int order);

double backwardRec(double x, double h, int order);

double centralRec(double x, double h, int order);

#endif