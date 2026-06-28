#include "../../include/pvi/eulerImplicito.h"

// Metodo de Euler Implicito
// Eq. S_{i+1} = S_i + dt * F(S_{i+1}, t_{i+1})
//
// O sistema implicito e resolvido pelo "resolutor" passado como parametro,
// que encapsula a solucao analitica (ou numerica) da equacao implicita
// para o PVI especifico.
std::vector<Estado> eulerImplicito(
    const ResolutorImplicito &resolver,
    const Estado             &S0,
    double                    t0,
    double                    dt,
    int                       nPassos,
    std::vector<double>      &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    estados.push_back(S0);
    tempos.push_back(t0);

    Estado Si = S0;
    double ti = t0;

    for (int i = 0; i < nPassos; i++)
    {
        // Chama o resolutor especifico do PVI para obter S_{i+1}
        Estado Snew = resolver(Si, ti, dt);

        ti += dt;
        estados.push_back(Snew);
        tempos.push_back(ti);

        Si = Snew;
    }

    return estados;
}