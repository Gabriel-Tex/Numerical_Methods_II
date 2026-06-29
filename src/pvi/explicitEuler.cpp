#include "../../include/pvi/explicitEuler.h"

std::vector<State> explicitEuler(
    const StateFunc &F,
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
        // S_{i+1} = S_i + dt * F(S_i, t_i)
        State Snew = addStates(Si, scaleState(dt, F(Si, ti)));
        ti += dt;
        states.push_back(Snew);
        times.push_back(ti);
        Si = Snew;
    }

    return states;
}