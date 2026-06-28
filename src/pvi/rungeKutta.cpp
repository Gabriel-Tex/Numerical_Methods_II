#include "../../include/pvi/rungeKutta.h"

// Runge-Kutta de 2a ordem
// P 1.1: S_bar_{i+1} = S_i + dt * F(S_i, t_i)            
// P 1.2: S_{i+1} = S_i + dt*(1/2*F1 + 1/2*F2)            
std::vector<Estado> rungeKutta2(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    estados.push_back(S0);
    tempos.push_back(t0);

    Estado Si = S0;
    double ti = t0;

    for (int i = 0; i < nPassos; i++)
    {
        // P 1.1: estimativa grosseira
        Estado F1 = F(Si, ti);
        Estado Sbar = somaEstado(Si, escalaEstado(dt, F1));

        // P 1.2: atualizacao melhorada
        Estado F2 = F(Sbar, ti + dt);
        Estado media = escalaEstado(0.5, somaEstado(F1, F2));
        Estado Snew = somaEstado(Si, escalaEstado(dt, media));

        ti += dt;
        estados.push_back(Snew);
        tempos.push_back(ti);

        Si = Snew;
    }

    return estados;
}

// Runge-Kutta de 3a ordem
// F1 = F(S_i, t_i)
// S_bar_{i+1/2} = S_i + (dt/2)*F1                        
// F2 = F(S_bar_{i+1/2}, t_i + dt/2)                       
// S_bar_{i+1}   = S_i + dt*(-F1 + 2*F2)                   
// F3 = F(S_bar_{i+1}, t_i + dt)                           
// S_{i+1} = S_i + dt*(1/6*F1 + 4/6*F2 + 1/6*F3)          
std::vector<Estado> rungeKutta3(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    estados.push_back(S0);
    tempos.push_back(t0);

    Estado Si = S0;
    double ti = t0;

    for (int i = 0; i < nPassos; i++)
    {
        // F1
        Estado F1 = F(Si, ti);

        // S_bar_{i+1/2} = S_i + (dt/2)*F1             
        Estado Sbar_meio = somaEstado(Si, escalaEstado(dt / 2.0, F1));

        // F2
        Estado F2 = F(Sbar_meio, ti + dt / 2.0); // [Eq. (49)]

        // S_bar_{i+1} = S_i + dt*(-F1 + 2*F2)          
        Estado mF1_2F2 = somaEstado(escalaEstado(-1.0, F1), escalaEstado(2.0, F2));
        Estado Sbar = somaEstado(Si, escalaEstado(dt, mF1_2F2));

        // F3
        Estado F3 = F(Sbar, ti + dt); // 

        // S_{i+1} = S_i + dt*(1/6*F1 + 4/6*F2 + 1/6*F3) 
        Estado media = somaEstado(
            somaEstado(escalaEstado(1.0 / 6.0, F1), escalaEstado(4.0 / 6.0, F2)),
            escalaEstado(1.0 / 6.0, F3));
        Estado Snew = somaEstado(Si, escalaEstado(dt, media));

        ti += dt;
        estados.push_back(Snew);
        tempos.push_back(ti);

        Si = Snew;
    }

    return estados;
}

// Runge-Kutta de 4a ordem 
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
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    estados.push_back(S0);
    tempos.push_back(t0);

    Estado Si = S0;
    double ti = t0;

    for (int i = 0; i < nPassos; i++)
    {
        // F1 = F(S_i, t_i)                               
        Estado F1 = F(Si, ti);

        // S2 = S_i + (dt/2)*F1                          
        Estado S2 = somaEstado(Si, escalaEstado(dt / 2.0, F1));

        // F2 = F(S2, t_i + dt/2)                         
        Estado F2 = F(S2, ti + dt / 2.0);

        // S3 = S_i + (dt/2)*F2                          
        Estado S3 = somaEstado(Si, escalaEstado(dt / 2.0, F2));

        // F3 = F(S3, t_i + dt/2)                         
        Estado F3 = F(S3, ti + dt / 2.0);

        // S4 = S_i + dt*F3                               
        Estado S4 = somaEstado(Si, escalaEstado(dt, F3));

        // F4 = F(S4, t_i + dt)                          
        Estado F4 = F(S4, ti + dt);

        // S_{i+1} = S_i + (dt/6)*(F1 + 2*F2 + 2*F3 + F4) 
        Estado soma = somaEstado(
            somaEstado(F1, escalaEstado(2.0, F2)),
            somaEstado(escalaEstado(2.0, F3), F4));
        Estado Snew = somaEstado(Si, escalaEstado(dt / 6.0, soma));

        ti += dt;
        estados.push_back(Snew);
        tempos.push_back(ti);

        Si = Snew;
    }

    return estados;
}