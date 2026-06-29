#include "../../include/pvi/implicitEuler.h"

std::vector<State> implicitEuler(
    const ImplicitSolver &solver,
    const State &S0,
    double t0,
    double dt,
    int nSteps,
    std::vector<double> &times)
{
    std::vector<State> states;
    times.clear();

    states.push_back(S0);
    times.push_back(t0);

    State Si = S0;
    double ti = t0;

    for (int i = 0; i < nSteps; i++)
    {
        State Snew = solver(Si, ti, dt);
        ti += dt;
        states.push_back(Snew);
        times.push_back(ti);
        Si = Snew;
    }

    return states;
}