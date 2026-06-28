#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "pviUtils.h"
#include <vector>

// Metodo de Runge-Kutta de 2a ordem
// P 1.1: S_bar_{i+1} = S_i + dt * F(S_i, t_i)                       
// P 1.2: S_{i+1} = S_i + dt*(1/2*F(S_i,t_i) + 1/2*F(S_bar_{i+1},t_{i+1}))  
std::vector<Estado> rungeKutta2(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos);

// Metodo de Runge-Kutta de 3a ordem 
// F1 = F(S_i, t_i)
// S_bar_{i+1/2} = S_i + (dt/2)*F1                                     
// F2 = F(S_bar_{i+1/2}, t_i + dt/2)                                   
// S_bar_{i+1} = S_i + dt*(-F1 + 2*F2)                                
// F3 = F(S_bar_{i+1}, t_i + dt)                                      
// S_{i+1} = S_i + dt*(1/6*F1 + 4/6*F2 + 1/6*F3)                    
std::vector<Estado> rungeKutta3(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos);

// Metodo de Runge-Kutta de 4a ordem
// F1 = F(S_i, t_i)                                                      
// S2 = S_i + (dt/2)*F1                                                 
// F2 = F(S2, t_i + dt/2)                                               
// S3 = S_i + (dt/2)*F2                                                  
// F3 = F(S3, t_i + dt/2)                                                
// S4 = S_i + dt*F3                                                      
// F4 = F(S4, t_i + dt)                                                 
// S_{i+1} = S_i + (dt/6)*(F1 + 2*F2 + 2*F3 + F4)                       
std::vector<Estado> rungeKutta4(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos);

#endif