#include "../../include/pvi/rungeKutta.h"

std::vector<State> rungeKutta2(
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
        State F1 = F(Si, ti);
        State Sbar = addStates(Si, scaleState(dt, F1));
        State F2 = F(Sbar, ti + dt);
        State Snew = addStates(Si, scaleState(dt, scaleState(0.5, addStates(F1, F2))));
        ti += dt;
        states.push_back(Snew);
        times.push_back(ti);
        Si = Snew;
    }

    return states;
}

std::vector<State> rungeKutta3(
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
        State F1 = F(Si, ti);
        State Sbar_mid = addStates(Si, scaleState(dt / 2.0, F1));
        State F2 = F(Sbar_mid, ti + dt / 2.0);
        State Sbar = addStates(Si, scaleState(dt, addStates(scaleState(-1.0, F1), scaleState(2.0, F2))));
        State F3 = F(Sbar, ti + dt);
        State avg = addStates(addStates(scaleState(1.0 / 6.0, F1), scaleState(4.0 / 6.0, F2)), scaleState(1.0 / 6.0, F3));
        State Snew = addStates(Si, scaleState(dt, avg));
        ti += dt;
        states.push_back(Snew);
        times.push_back(ti);
        Si = Snew;
    }

    return states;
}

std::vector<State> rungeKutta4(
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
        State F1 = F(Si, ti);
        State S2 = addStates(Si, scaleState(dt / 2.0, F1));
        State F2 = F(S2, ti + dt / 2.0);
        State S3 = addStates(Si, scaleState(dt / 2.0, F2));
        State F3 = F(S3, ti + dt / 2.0);
        State S4 = addStates(Si, scaleState(dt, F3));
        State F4 = F(S4, ti + dt);
        State sum = addStates(addStates(F1, scaleState(2.0, F2)), addStates(scaleState(2.0, F3), F4));
        State Snew = addStates(Si, scaleState(dt / 6.0, sum));
        ti += dt;
        states.push_back(Snew);
        times.push_back(ti);
        Si = Snew;
    }

    return states;
}