#include "../../include/pvi/eulerExplicito.h"

//  Metodo de Euler Explicito
// Eq. (12): S_{i+1} = S_i + Delta_t * F(S_i, t_i)
std::vector<Estado> eulerExplicito(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    // Condicao inicial
    estados.push_back(S0);
    tempos.push_back(t0);

    Estado Si = S0;
    double ti = t0;

    for (int i = 0; i < nPassos; i++)
    {
        // Eq. (12): S_{i+1} = S_i + dt * F(S_i, t_i)
        Estado dS = F(Si, ti);
        Estado Snew = somaEstado(Si, escalaEstado(dt, dS));

        ti += dt;
        estados.push_back(Snew);
        tempos.push_back(ti);

        Si = Snew;
    }

    return estados;
}