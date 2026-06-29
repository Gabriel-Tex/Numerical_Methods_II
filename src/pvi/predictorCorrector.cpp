#include "../../include/pvi/predictorCorrector.h"
#include "../../include/pvi/rungeKutta.h"
#include <cmath>

std::vector<State> adamsBashforth2(
    const StateFunc &F, const State &S0, double t0, double dt,
    int nSteps, double eps, int maxIter, std::vector<double> &times)
{
    std::vector<State> states;
    times.clear();

    std::vector<double> tInit;
    auto init = rungeKutta2(F, S0, t0, dt, 1, tInit);
    states.push_back(init[0]);
    times.push_back(tInit[0]);
    states.push_back(init[1]);
    times.push_back(tInit[1]);

    std::vector<State> deriv;
    deriv.push_back(F(states[0], times[0]));
    deriv.push_back(F(states[1], times[1]));

    for (int i = 1; i < nSteps; i++)
    {
        const State &Si = states[i];
        double ti = times[i];
        const State &Fi = deriv[i];
        const State &Fi1 = deriv[i - 1];

        // Predição
        State Sbar = addStates(Si, scaleState(dt / 2.0, addStates(scaleState(-1.0, Fi1), scaleState(3.0, Fi))));

        // Correção iterativa
        State Snew = Sbar;
        for (int it = 0; it < maxIter; it++)
        {
            State Sold = Snew;
            State Fbar = F(Sbar, ti + dt);
            Snew = addStates(Si, scaleState(dt / 2.0, addStates(Fi, Fbar)));
            if (relativeError(Snew, Sold) < eps)
                break;
            Sbar = Snew;
        }

        times.push_back(ti + dt);
        states.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return states;
}

std::vector<State> adamsBashforthMoulton3(
    const StateFunc &F, const State &S0, double t0, double dt,
    int nSteps, double eps, int maxIter, std::vector<double> &times)
{
    std::vector<State> states;
    times.clear();

    std::vector<double> tInit;
    auto init = rungeKutta3(F, S0, t0, dt, 2, tInit);
    for (int k = 0; k <= 2; k++)
    {
        states.push_back(init[k]);
        times.push_back(tInit[k]);
    }

    std::vector<State> deriv;
    for (int k = 0; k <= 2; k++)
        deriv.push_back(F(states[k], times[k]));

    for (int i = 2; i < nSteps; i++)
    {
        const State &Si = states[i];
        double ti = times[i];
        const State &Fi = deriv[i], &Fi1 = deriv[i - 1], &Fi2 = deriv[i - 2];

        // Predição
        State Sbar = addStates(Si, scaleState(dt / 12.0,
                                              addStates(addStates(scaleState(5.0, Fi2), scaleState(-16.0, Fi1)), scaleState(23.0, Fi))));

        // Correção iterativa
        State Snew = Sbar;
        for (int it = 0; it < maxIter; it++)
        {
            State Sold = Snew;
            State Fbar = F(Sbar, ti + dt);
            Snew = addStates(Si, scaleState(dt / 12.0,
                                            addStates(addStates(scaleState(-1.0, Fi1), scaleState(8.0, Fi)), scaleState(5.0, Fbar))));
            if (relativeError(Snew, Sold) < eps)
                break;
            Sbar = Snew;
        }

        times.push_back(ti + dt);
        states.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return states;
}

std::vector<State> adamsBashforthMoulton4(
    const StateFunc &F, const State &S0, double t0, double dt,
    int nSteps, double eps, int maxIter, std::vector<double> &times)
{
    std::vector<State> states;
    times.clear();

    std::vector<double> tInit;
    auto init = rungeKutta4(F, S0, t0, dt, 3, tInit);
    for (int k = 0; k <= 3; k++)
    {
        states.push_back(init[k]);
        times.push_back(tInit[k]);
    }

    std::vector<State> deriv;
    for (int k = 0; k <= 3; k++)
        deriv.push_back(F(states[k], times[k]));

    for (int i = 3; i < nSteps; i++)
    {
        const State &Si = states[i];
        double ti = times[i];
        const State &Fi = deriv[i], &Fi1 = deriv[i - 1], &Fi2 = deriv[i - 2], &Fi3 = deriv[i - 3];

        // Predição
        State Sbar = addStates(Si, scaleState(dt / 24.0,
                                              addStates(addStates(scaleState(-9.0, Fi3), scaleState(37.0, Fi2)),
                                                        addStates(scaleState(-59.0, Fi1), scaleState(55.0, Fi)))));

        // Correção iterativa
        State Snew = Sbar;
        for (int it = 0; it < maxIter; it++)
        {
            State Sold = Snew;
            State Fbar = F(Sbar, ti + dt);
            Snew = addStates(Si, scaleState(dt / 24.0,
                                            addStates(addStates(scaleState(1.0, Fi2), scaleState(-5.0, Fi1)),
                                                      addStates(scaleState(19.0, Fi), scaleState(9.0, Fbar)))));
            if (relativeError(Snew, Sold) < eps)
                break;
            Sbar = Snew;
        }

        times.push_back(ti + dt);
        states.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return states;
}