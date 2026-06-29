#include "../../include/pvi/menuPVI.h"
#include "../../include/pvi/pviUtils.h"
#include "../../include/pvi/explicitEuler.h"
#include "../../include/pvi/implicitEuler.h"
#include "../../include/pvi/rungeKutta.h"
#include "../../include/pvi/predictorCorrector.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// ---- utilitários locais ----

static void pause()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int readInt(const std::string &prompt, int min, int max)
{
    int v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v && v >= min && v <= max)
            return v;
        std::cout << "  Opcao invalida. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static double readDouble(const std::string &prompt)
{
    double v;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> v)
            return v;
        std::cout << "  Valor invalido. Tente novamente.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// PVI-1: dy/dt = (2/3)*y
static State F_IVP1(const State &S, double)
{
    return {(2.0 / 3.0) * S[0]};
}

// PVI-2: dv/dt = -g - (k/m)*v,  dy/dt = v
static const double G_IVP2 = 10.0;
static const double K_IVP2 = 0.25;
static const double M_IVP2 = 2.0;

static State F_IVP2(const State &S, double)
{
    return {-G_IVP2 - (K_IVP2 / M_IVP2) * S[0], S[0]};
}

// Resolutor implicito para o PVI-2
static State implicitSolverIVP2(const State &Si, double, double dt)
{
    double den = M_IVP2 + K_IVP2 * dt;
    double vNew = (M_IVP2 / den) * (Si[0] - G_IVP2 * dt);
    return {vNew, Si[1] + dt * vNew};
}

// Solucao exata do PVI-2
static double exactV_IVP2(double t, double t0, double v0)
{
    double gm_k = G_IVP2 * M_IVP2 / K_IVP2;
    return -gm_k + (v0 + gm_k) * std::exp(-(K_IVP2 / M_IVP2) * (t - t0));
}

static double exactY_IVP2(double t, double t0, double v0, double y0)
{
    double gm_k = G_IVP2 * M_IVP2 / K_IVP2;
    double mk = M_IVP2 / K_IVP2;
    return y0 - gm_k * (t - t0) - (v0 + gm_k) * mk * (std::exp(-(K_IVP2 / M_IVP2) * (t - t0)) - 1.0);
}

// Exibe tabela v, y, solucao exata e erro para PVI-2
static void printIVP2Table(
    const std::vector<State> &states,
    const std::vector<double> &times,
    double t0, double v0, double y0)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n"
              << std::setw(8) << "t"
              << std::setw(14) << "v_aprox"
              << std::setw(14) << "y_aprox"
              << std::setw(14) << "v_exato"
              << std::setw(14) << "y_exato"
              << std::setw(12) << "err_v(%)"
              << std::setw(12) << "err_y(%)\n";
    std::cout << std::string(88, '-') << "\n";

    for (int i = 0; i < (int)times.size(); i++)
    {
        double t = times[i];
        double vAp = states[i][0], yAp = states[i][1];
        double vEx = exactV_IVP2(t, t0, v0);
        double yEx = exactY_IVP2(t, t0, v0, y0);
        double ev = (std::fabs(vEx) > 1e-15) ? (vAp - vEx) / vEx * 100.0 : 0.0;
        double ey = (std::fabs(yEx) > 1e-15) ? (yAp - yEx) / yEx * 100.0 : 0.0;
        std::cout << std::setw(8) << t << std::setw(14) << vAp << std::setw(14) << yAp
                  << std::setw(14) << vEx << std::setw(14) << yEx
                  << std::setw(12) << ev << std::setw(12) << ey << "\n";
    }
}

// Analisa trajetoria (altura máxima, tempo de impacto, velocidade de impacto)
static void analyzeTrajectory(
    const std::vector<State> &states,
    const std::vector<double> &times,
    double &yMax, double &tYmax, double &tTotal, double &vImpact)
{
    yMax = states[0][1];
    tYmax = times[0];
    tTotal = times.back();
    vImpact = states.back()[0];

    for (int i = 0; i < (int)times.size(); i++)
    {
        if (states[i][1] > yMax)
        {
            yMax = states[i][1];
            tYmax = times[i];
        }
        if (i > 0 && states[i][1] <= 0.0 && states[i - 1][1] > 0.0)
        {
            double y0 = states[i - 1][1], y1 = states[i][1];
            double t0l = times[i - 1], t1 = times[i];
            tTotal = t0l + (0.0 - y0) / (y1 - y0) * (t1 - t0l);
            vImpact = states[i][0];
        }
    }
}

static const double T0_TASK = 0.0;
static const double V0_TASK = 5.0;
static const double Y0_TASK = 200.0;

static void runIVP2WithMultipleDt(
    const std::string &methodName,
    std::function<std::vector<State>(double, int, std::vector<double> &)> runMethod)
{
    std::vector<double> dts = {0.1, 0.01, 0.001, 0.0001};
    double tFinal = 100.0;

    for (double dt : dts)
    {
        int nSteps = (int)(tFinal / dt);
        std::vector<double> times;
        auto states = runMethod(dt, nSteps, times);

        int iImpact = (int)times.size() - 1;
        for (int i = 1; i < (int)states.size(); i++)
            if (states[i][1] <= 0.0)
            {
                iImpact = i;
                break;
            }

        std::vector<State> stC(states.begin(), states.begin() + iImpact + 1);
        std::vector<double> tC(times.begin(), times.begin() + iImpact + 1);

        double yMax, tYmax, tTotal, vImpact;
        analyzeTrajectory(stC, tC, yMax, tYmax, tTotal, vImpact);

        std::cout << "\n--- " << methodName << "  (dt = " << dt << ") ---\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  a) Altura maxima:              y_max  = " << yMax << " m\n";
        std::cout << "  b) Tempo ate altura maxima:    t_ymax = " << tYmax << " s\n";
        std::cout << "  c) Tempo total ate o mar:      t_tot  = " << tTotal << " s\n";
        std::cout << "  d) Velocidade no impacto:      v_imp  = " << vImpact << " m/s\n";
    }
}

// ---- Sub-menus das aulas ----

static void menuEuler()
{
    std::cout << "   Euler Explicito / Implicito            \n";
    std::cout << "  1 - PVI-1 (demo passo a passo)         \n";
    std::cout << "  2 - Tarefa 16: PVI-2 com varios dt     \n";
    std::cout << "  0 - Voltar                             \n";

    int op = readInt("Opcao: ", 0, 2);

    if (op == 1)
    {
        double dt = readDouble("Delta t: ");
        int nSteps = readInt("Numero de passos: ", 1, 100000);

        std::cout << "\n--- Euler Explicito ---\n";
        std::vector<double> times;
        auto states = explicitEuler(F_IVP1, {2.0}, 0.0, dt, nSteps, times);
        std::cout << std::fixed << std::setprecision(6);
        for (int i = 0; i < (int)times.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * times[i]);
            double err = std::fabs((states[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << times[i] << "  y_aprox=" << states[i][0]
                      << "  y_exato=" << yEx << "  erro=" << err << "%\n";
        }

        std::cout << "\n--- Euler Implicito ---\n";
        // Resolutor analitico para PVI-1: y_{i+1} = y_i / (1 - 2/3 * dt)
        ImplicitSolver solvePVI1 = [](const State &Si, double, double dt) -> State
        {
            return {Si[0] / (1.0 - (2.0 / 3.0) * dt)};
        };
        std::vector<double> times2;
        auto states2 = implicitEuler(solvePVI1, {2.0}, 0.0, dt, nSteps, times2);
        for (int i = 0; i < (int)times2.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * times2[i]);
            double err = std::fabs((states2[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << times2[i] << "  y_aprox=" << states2[i][0]
                      << "  y_exato=" << yEx << "  erro=" << err << "%\n";
        }
        pause();
    }
    else if (op == 2)
    {
        std::cout << "\n=== Tarefa 16: Euler Explicito e Implicito ===\n";
        std::cout << "PVI-2: t0=0, v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";

        State S0 = {V0_TASK, Y0_TASK};

        std::cout << "\n*** Euler Explicito ***\n";
        runIVP2WithMultipleDt("Euler Explicito",
                              [&](double dt, int n, std::vector<double> &t)
                              {
                                  return explicitEuler(F_IVP2, S0, T0_TASK, dt, n, t);
                              });

        std::cout << "\n*** Euler Implicito ***\n";
        runIVP2WithMultipleDt("Euler Implicito",
                              [&](double dt, int n, std::vector<double> &t)
                              {
                                  return implicitEuler(implicitSolverIVP2, S0, T0_TASK, dt, n, t);
                              });

        pause();
    }
}

static void menuRungeKutta()
{
    std::cout << "   Runge-Kutta                            \n";
    std::cout << "  1 - RK2 ou RK3 no PVI-1 (demo)          \n";
    std::cout << "  2 - Tarefa 17: PVI-2 com RK3            \n";
    std::cout << "  0 - Voltar                              \n";

    int op = readInt("Opcao: ", 0, 2);

    if (op == 1)
    {
        int order = readInt("Ordem (2 ou 3): ", 2, 3);
        double dt = readDouble("Delta t: ");
        int nSteps = readInt("Numero de passos: ", 1, 100000);

        std::vector<double> times;
        std::vector<State> states;
        if (order == 2)
            states = rungeKutta2(F_IVP1, {2.0}, 0.0, dt, nSteps, times);
        else
            states = rungeKutta3(F_IVP1, {2.0}, 0.0, dt, nSteps, times);

        std::cout << std::fixed << std::setprecision(6);
        for (int i = 0; i < (int)times.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * times[i]);
            double err = std::fabs((states[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << times[i] << "  y=" << states[i][0]
                      << "  y_ex=" << yEx << "  erro=" << err << "%\n";
        }
        pause();
    }
    else if (op == 2)
    {
        std::cout << "\n=== Tarefa 17: Runge-Kutta 3a Ordem ===\n";
        State S0 = {V0_TASK, Y0_TASK};
        runIVP2WithMultipleDt("RK3",
                              [&](double dt, int n, std::vector<double> &t)
                              {
                                  return rungeKutta3(F_IVP2, S0, T0_TASK, dt, n, t);
                              });
        pause();
    }
}

static void menuPredictorCorrector()
{
    std::cout << "   Preditor-Corretor                      \n";
    std::cout << "  1 - AB2  no PVI-1 (demo)               \n";
    std::cout << "  2 - ABM3 no PVI-1 (demo)               \n";
    std::cout << "  3 - ABM4 no PVI-1 (demo)               \n";
    std::cout << "  4 - Tarefa 18: ABM4 no PVI-2           \n";
    std::cout << "  0 - Voltar                             \n";

    int op = readInt("Opcao: ", 0, 4);

    if (op >= 1 && op <= 3)
    {
        double dt = readDouble("Delta t: ");
        int nSteps = readInt("Numero de passos: ", 1, 100000);
        double eps = 1e-10;
        int maxIter = 50;

        std::vector<double> times;
        std::vector<State> states;
        std::string name;

        if (op == 1)
        {
            states = adamsBashforth2(F_IVP1, {2.0}, 0.0, dt, nSteps, eps, maxIter, times);
            name = "AB2";
        }
        else if (op == 2)
        {
            states = adamsBashforthMoulton3(F_IVP1, {2.0}, 0.0, dt, nSteps, eps, maxIter, times);
            name = "ABM3";
        }
        else
        {
            states = adamsBashforthMoulton4(F_IVP1, {2.0}, 0.0, dt, nSteps, eps, maxIter, times);
            name = "ABM4";
        }

        std::cout << "\n--- " << name << " no PVI-1 ---\n";
        std::cout << std::fixed << std::setprecision(8);
        for (int i = 0; i < (int)times.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * times[i]);
            double err = std::fabs((states[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << times[i] << "  y=" << states[i][0]
                      << "  y_ex=" << yEx << "  erro=" << err << "%\n";
        }
        pause();
    }
    else if (op == 4)
    {
        std::cout << "\n=== Tarefa 18: ABM 4a ordem no PVI-2 ===\n";
        State S0 = {V0_TASK, Y0_TASK};
        double eps = 1e-10;
        int maxIter = 50;

        runIVP2WithMultipleDt("ABM4",
                              [&](double dt, int n, std::vector<double> &t)
                              {
                                  return adamsBashforthMoulton4(F_IVP2, S0, T0_TASK, dt, n, maxIter, eps, t);
                              });
        pause();
    }
}

void menuPVI()
{
    while (true)
    {
        std::cout << "   PVI - Problemas de Valor Inicial       \n";
        std::cout << "  1 - Euler Explicito / Implicito         \n";
        std::cout << "  2 - Runge-Kutta 2a e 3a ord.            \n";
        std::cout << "  3 - Preditor-Corretor                   \n";
        std::cout << "  0 - Voltar                              \n";

        int op = readInt("Opcao: ", 0, 3);
        if (op == 0)
            return;
        if (op == 1)
            menuEuler();
        if (op == 2)
            menuRungeKutta();
        if (op == 3)
            menuPredictorCorrector();
    }
}