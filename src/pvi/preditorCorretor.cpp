#include "../../include/pvi/preditorCorretor.h"
#include "../../include/pvi/rungeKutta.h"
#include <cmath>

// Adams-Bashforth de 2a ordem (Predicao-Correcao)
//
// Inicializacao: S1 por RK2
// Predicao  [Eq. 19]: S_bar = S_i + (dt/2)*(-F_{i-1} + 3*F_i)
// Correcao  [Eq. 20]: S_{i+1} = S_i + (dt/2)*(F_i + F(S_bar, t_{i+1}))
std::vector<Estado> adamsBashforth2(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    // --- Fase de inicializacao: S1 por RK2 ---
    std::vector<double> tInit;
    auto initStates = rungeKutta2(F, S0, t0, dt, 1, tInit);

    estados.push_back(initStates[0]);
    tempos.push_back(tInit[0]);
    estados.push_back(initStates[1]);
    tempos.push_back(tInit[1]);

    // Pre-calcula as derivadas nos passos inicializados
    std::vector<Estado> deriv;
    deriv.push_back(F(estados[0], tempos[0])); // F_{i-1} = F_0
    deriv.push_back(F(estados[1], tempos[1])); // F_i     = F_1

    // --- Fase de predicao e correcao ---
    for (int i = 1; i < nPassos; i++)
    {
        const Estado &Si = estados[i];
        double ti = tempos[i];
        const Estado &Fi = deriv[i];
        const Estado &Fi_1 = deriv[i - 1];

        // Predicao: S_bar = S_i + (dt/2)*(-F_{i-1} + 3*F_i)
        Estado passo_pred = somaEstado(escalaEstado(-1.0, Fi_1), escalaEstado(3.0, Fi));
        Estado Sbar = somaEstado(Si, escalaEstado(dt / 2.0, passo_pred));

        // Correcao iterativa 
        Estado Snew = Sbar;
        for (int iter = 0; iter < maxIter; iter++)
        {
            Estado Snew_ant = Snew;

            Estado Fbar = F(Sbar, ti + dt);
            Estado media = escalaEstado(0.5, somaEstado(Fi, Fbar));
            Snew = somaEstado(Si, escalaEstado(dt, media));

            if (erroRelativo(Snew, Snew_ant) < eps)
                break;

            Sbar = Snew;
        }

        tempos.push_back(ti + dt);
        estados.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return estados;
}

// Adams-Bashforth-Moulton de 3a ordem
//
// Inicializacao: S1, S2 por RK3
// Predicao  [Eq. 43]: S_bar = S_i + (dt/12)*(5*F_{i-2} - 16*F_{i-1} + 23*F_i)
// Correcao  [Eq. 44]: S_{i+1} = S_i + (dt/12)*(-F_{i-1} + 8*F_i + 5*F(S_bar,t_{i+1}))
std::vector<Estado> adamsBashforthMoulton3(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    // --- Fase de inicializacao: S1 e S2 por RK3 ---
    std::vector<double> tInit;
    auto initStates = rungeKutta3(F, S0, t0, dt, 2, tInit);

    for (int k = 0; k <= 2; k++)
    {
        estados.push_back(initStates[k]);
        tempos.push_back(tInit[k]);
    }

    // Pre-calcula derivadas
    std::vector<Estado> deriv;
    for (int k = 0; k <= 2; k++)
        deriv.push_back(F(estados[k], tempos[k]));

    // --- Fase de predicao e correcao (comeca em i = 2) ---
    for (int i = 2; i < nPassos; i++)
    {
        const Estado &Si = estados[i];
        double ti = tempos[i];
        const Estado &Fi = deriv[i];
        const Estado &Fi_1 = deriv[i - 1];
        const Estado &Fi_2 = deriv[i - 2];

        // Predicao:
        // S_bar = S_i + (dt/12)*(5*F_{i-2} - 16*F_{i-1} + 23*F_i)
        Estado passo_pred = somaEstado(
            somaEstado(escalaEstado(5.0, Fi_2), escalaEstado(-16.0, Fi_1)),
            escalaEstado(23.0, Fi));
        Estado Sbar = somaEstado(Si, escalaEstado(dt / 12.0, passo_pred));

        // Correcao iterativa 
        Estado Snew = Sbar;
        for (int iter = 0; iter < maxIter; iter++)
        {
            Estado Snew_ant = Snew;

            Estado Fbar = F(Sbar, ti + dt);
            // (-F_{i-1} + 8*F_i + 5*F_bar)
            Estado passo_cor = somaEstado(
                somaEstado(escalaEstado(-1.0, Fi_1), escalaEstado(8.0, Fi)),
                escalaEstado(5.0, Fbar));
            Snew = somaEstado(Si, escalaEstado(dt / 12.0, passo_cor));

            if (erroRelativo(Snew, Snew_ant) < eps)
                break;

            Sbar = Snew;
        }

        tempos.push_back(ti + dt);
        estados.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return estados;
}

// Adams-Bashforth-Moulton de 4a ordem
//
// Inicializacao: S1, S2, S3 por RK4
// Predicao:
//   S_bar = S_i + (dt/24)*(-9*F_{i-3} + 37*F_{i-2} - 59*F_{i-1} + 55*F_i)
// Correcao:
//   S_{i+1} = S_i + (dt/24)*(F_{i-2} - 5*F_{i-1} + 19*F_i + 9*F_bar)
std::vector<Estado> adamsBashforthMoulton4(
    const FuncaoF &F,
    const Estado &S0,
    double t0,
    double dt,
    int nPassos,
    double eps,
    int maxIter,
    std::vector<double> &tempos)
{
    std::vector<Estado> estados;
    tempos.clear();

    // --- Fase de inicializacao: S1, S2, S3 por RK4 ---
    std::vector<double> tInit;
    auto initStates = rungeKutta4(F, S0, t0, dt, 3, tInit);

    for (int k = 0; k <= 3; k++)
    {
        estados.push_back(initStates[k]);
        tempos.push_back(tInit[k]);
    }

    // Pre-calcula derivadas
    std::vector<Estado> deriv;
    for (int k = 0; k <= 3; k++)
        deriv.push_back(F(estados[k], tempos[k]));

    // --- Fase de predicao e correcao (comeca em i = 3) ---
    for (int i = 3; i < nPassos; i++)
    {
        const Estado &Si = estados[i];
        double ti = tempos[i];
        const Estado &Fi = deriv[i];
        const Estado &Fi_1 = deriv[i - 1];
        const Estado &Fi_2 = deriv[i - 2];
        const Estado &Fi_3 = deriv[i - 3];

        // Predicao: S_bar = S_i + (dt/24)*(-9*F_{i-3}+37*F_{i-2}-59*F_{i-1}+55*F_i)
        Estado passo_pred = somaEstado(
            somaEstado(escalaEstado(-9.0, Fi_3), escalaEstado(37.0, Fi_2)),
            somaEstado(escalaEstado(-59.0, Fi_1), escalaEstado(55.0, Fi)));
        Estado Sbar = somaEstado(Si, escalaEstado(dt / 24.0, passo_pred));

        // Correcao iterativa
        Estado Snew = Sbar;
        for (int iter = 0; iter < maxIter; iter++)
        {
            Estado Snew_ant = Snew;

            Estado Fbar = F(Sbar, ti + dt);
            // (F_{i-2} - 5*F_{i-1} + 19*F_i + 9*F_bar)
            Estado passo_cor = somaEstado(
                somaEstado(escalaEstado(1.0, Fi_2), escalaEstado(-5.0, Fi_1)),
                somaEstado(escalaEstado(19.0, Fi), escalaEstado(9.0, Fbar)));
            Snew = somaEstado(Si, escalaEstado(dt / 24.0, passo_cor));

            if (erroRelativo(Snew, Snew_ant) < eps)
                break;

            Sbar = Snew;
        }

        tempos.push_back(ti + dt);
        estados.push_back(Snew);
        deriv.push_back(F(Snew, ti + dt));
    }

    return estados;
}