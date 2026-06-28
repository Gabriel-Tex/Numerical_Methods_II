#include "../../include/pvi/menuPVI.h"
#include "../../include/pvi/pviUtils.h"
#include "../../include/pvi/eulerExplicito.h"
#include "../../include/pvi/eulerImplicito.h"
#include "../../include/pvi/rungeKutta.h"
#include "../../include/pvi/preditorCorretor.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// ----- Utilitarios locais -------

static void pausarPVI()
{
    std::cout << "\n[Enter para continuar...]";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

static int lerIntPVI(const std::string &prompt, int min, int max)
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

static double lerDoublePVI(const std::string &prompt)
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

// PVI-1: dy/dt = (2/3)*y,  solucao exata: y(t) = y0 * exp(2/3*(t-t0))
static Estado F_PVI1(const Estado &S, double /*t*/)
{
    return {(2.0 / 3.0) * S[0]};
}

// PVI-2 (particula em queda com resistencia do ar):
//   dv/dt = -g - (k/m)*v
//   dy/dt = v
// Estado S = {v, y},  parametros: g=10, k=0.25, m=2
static const double G_PVI2 = 10.0;
static const double K_PVI2 = 0.25;
static const double M_PVI2 = 2.0;

static Estado F_PVI2(const Estado &S, double /*t*/)
{
    double v = S[0];
    double dvdt = -G_PVI2 - (K_PVI2 / M_PVI2) * v;
    double dydt = v;
    return {dvdt, dydt};
}

// Resolutor analitico do Euler Implicito para o PVI-2
// Sistema linear: [1 + (k/m)*dt  0 ] [v_{i+1}]   [v_i - g*dt]
//                 [-dt            1 ] [y_{i+1}] = [y_i       ]
// Solucao:
//   v_{i+1} = m/(m + k*dt) * (v_i - g*dt)
//   y_{i+1} = y_i + dt/(m + k*dt)*m * (v_i - g*dt)
//           = y_i + dt * v_{i+1}
static Estado resolutorImplicitoPVI2(const Estado &Si, double /*ti*/, double dt)
{
    double v_i = Si[0];
    double y_i = Si[1];
    double den = M_PVI2 + K_PVI2 * dt;
    double v_new = (M_PVI2 / den) * (v_i - G_PVI2 * dt);
    double y_new = y_i + dt * v_new;
    return {v_new, y_new};
}

// Solucao exata do PVI-2:
//   v(t) = -g*m/k + (v0 + g*m/k)*exp(-k/m*(t-t0))
//   y(t) = y0 - (g*m/k)*(t-t0) - (v0 + g*m/k)*(m/k)*(exp(-k/m*(t-t0)) - 1)
static double vExatoPVI2(double t, double t0, double v0)
{
    double gm_k = G_PVI2 * M_PVI2 / K_PVI2;
    return -gm_k + (v0 + gm_k) * std::exp(-(K_PVI2 / M_PVI2) * (t - t0));
}
static double yExatoPVI2(double t, double t0, double v0, double y0)
{
    double gm_k = G_PVI2 * M_PVI2 / K_PVI2;
    double mk = M_PVI2 / K_PVI2;
    return y0 - gm_k * (t - t0) - (v0 + gm_k) * mk * (std::exp(-(K_PVI2 / M_PVI2) * (t - t0)) - 1.0);
}

// Funcoes auxiliares de exibicao

// Exibe tabela com v, y, solucao exata e erro (para PVI-2)
static void exibeTabelaPVI2(
    const std::vector<Estado> &estados,
    const std::vector<double> &tempos,
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

    for (int i = 0; i < (int)tempos.size(); i++)
    {
        double t = tempos[i];
        double v_ap = estados[i][0];
        double y_ap = estados[i][1];
        double v_ex = vExatoPVI2(t, t0, v0);
        double y_ex = yExatoPVI2(t, t0, v0, y0);
        double ev = (std::fabs(v_ex) > 1e-15) ? (v_ap - v_ex) / v_ex * 100.0 : 0.0;
        double ey = (std::fabs(y_ex) > 1e-15) ? (y_ap - y_ex) / y_ex * 100.0 : 0.0;

        std::cout << std::setw(8) << t
                  << std::setw(14) << v_ap
                  << std::setw(14) << y_ap
                  << std::setw(14) << v_ex
                  << std::setw(14) << y_ex
                  << std::setw(12) << ev
                  << std::setw(12) << ey << "\n";
    }
}

// Extrai y_max, t_ymax, t_total e v_impacto de uma solucao do PVI-2
static void analisaTrajetoria(
    const std::vector<Estado> &estados,
    const std::vector<double> &tempos,
    double &yMax, double &tYmax, double &tTotal, double &vImpacto)
{
    yMax = estados[0][1];
    tYmax = tempos[0];
    tTotal = tempos.back();
    vImpacto = estados.back()[0];

    for (int i = 0; i < (int)tempos.size(); i++)
    {
        double y = estados[i][1];
        if (y > yMax)
        {
            yMax = y;
            tYmax = tempos[i];
        }
        // Detecta cruzamento y = 0
        if (i > 0 && estados[i][1] <= 0.0 && estados[i - 1][1] > 0.0)
        {
            // Interpolacao linear para achar t_total exato
            double y0 = estados[i - 1][1], y1 = estados[i][1];
            double t0 = tempos[i - 1], t1 = tempos[i];
            tTotal = t0 + (0.0 - y0) / (y1 - y0) * (t1 - t0);
            vImpacto = estados[i][0];
        }
    }
}

// Condicoes iniciais das Tarefas 16, 17, 18: t0=0, v0=5, y0=200, k=0.25, m=2
static const double T0_TAREFA = 0.0;
static const double V0_TAREFA = 5.0;
static const double Y0_TAREFA = 200.0;

static void executaTarefaPVI2(const std::string &nomeMetodo, double dt,
                              const std::vector<Estado> &estados,
                              const std::vector<double> &tempos)
{
    double yMax, tYmax, tTotal, vImpacto;
    analisaTrajetoria(estados, tempos, yMax, tYmax, tTotal, vImpacto);

    std::cout << "\n--- " << nomeMetodo << "  (dt = " << dt << ") ---\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  a) Altura maxima:              y_max  = " << yMax << " m\n";
    std::cout << "  b) Tempo ate altura maxima:    t_ymax = " << tYmax << " s\n";
    std::cout << "  c) Tempo total ate o mar:      t_tot  = " << tTotal << " s\n";
    std::cout << "  d) Velocidade no impacto:      v_imp  = " << vImpacto << " m/s\n";
}

// Roda o metodo escolhido para varios valores de dt e exibe os resultados
static void rodaPVI2ComVariosDt(
    const std::string &nomeMetodo,
    std::function<std::vector<Estado>(double dt, int nPassos, std::vector<double> &)> executaMetodo)
{
    // Tabela de dt pedida nas tarefas
    std::vector<double> dts = {0.1, 0.01, 0.001, 0.0001};
    double tFinal = 100.0; // tempo suficientemente longo para a queda

    for (double dt : dts)
    {
        int nPassos = (int)(tFinal / dt);
        std::vector<double> tempos;
        auto estados = executaMetodo(dt, nPassos, tempos);

        // Encontra o passo do impacto (y <= 0) para cortar a saida
        int iImpacto = (int)tempos.size() - 1;
        for (int i = 1; i < (int)estados.size(); i++)
            if (estados[i][1] <= 0.0)
            {
                iImpacto = i;
                break;
            }

        std::vector<Estado> estCort(estados.begin(), estados.begin() + iImpacto + 1);
        std::vector<double> tCort(tempos.begin(), tempos.begin() + iImpacto + 1);

        executaTarefaPVI2(nomeMetodo, dt, estCort, tCort);
    }
}

// ----- Euler Explicito e Implicito -----
static void menuAula24()
{
    std::cout << "\n══════════════════════════════════════\n";
    std::cout << "   Euler Explicito/Implicito \n";
    std::cout << "══════════════════════════════════════\n";
    std::cout << "  1 - PVI-1 (demo passo a passo)\n";
    std::cout << "  2 - PVI-2 com varios dt\n";
    std::cout << "  0 - Voltar\n";

    int op = lerIntPVI("Opcao: ", 0, 2);

    if (op == 1)
    {
        // PVI-1: y(t) = 2*exp(2/3*t), demo com dt = 0.5
        double dt = lerDoublePVI("Delta t: ");
        int nPassos = lerIntPVI("Numero de passos: ", 1, 100000);

        std::cout << "\n--- Euler Explicito ---\n";
        std::vector<double> tempos;
        auto est = eulerExplicito(F_PVI1, {2.0}, 0.0, dt, nPassos, tempos);
        std::cout << std::fixed << std::setprecision(6);
        for (int i = 0; i < (int)tempos.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * tempos[i]);
            double err = std::fabs((est[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << tempos[i]
                      << "  y_aprox=" << est[i][0]
                      << "  y_exato=" << yEx
                      << "  erro=" << err << "%\n";
        }

        std::cout << "\n--- Euler Implicito ---\n";
        // Resolutor analitico para PVI-1: y_{i+1} = y_i / (1 - 2/3*dt)
        ResolutorImplicito resolvePVI1 = [](const Estado &Si, double /*ti*/, double dt) -> Estado
        {
            return {Si[0] / (1.0 - (2.0 / 3.0) * dt)};
        };
        std::vector<double> tempos2;
        auto est2 = eulerImplicito(resolvePVI1, {2.0}, 0.0, dt, nPassos, tempos2);
        for (int i = 0; i < (int)tempos2.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * tempos2[i]);
            double err = std::fabs((est2[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << tempos2[i]
                      << "  y_aprox=" << est2[i][0]
                      << "  y_exato=" << yEx
                      << "  erro=" << err << "%\n";
        }
        pausarPVI();
    }
    else if (op == 2)
    {
        // PVI-2 com t0=0, v0=5, y0=200, k=0.25, m=2
        std::cout << "\n=== Tarefa 16: Euler Explicito e Implicito ===\n";
        std::cout << "PVI-2: t0=0, v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";

        Estado S0 = {V0_TAREFA, Y0_TAREFA};

        std::cout << "\n*** Euler Explicito ***\n";
        rodaPVI2ComVariosDt("Euler Explicito",
                            [&](double dt, int nPassos, std::vector<double> &t)
                            {
                                return eulerExplicito(F_PVI2, S0, T0_TAREFA, dt, nPassos, t);
                            });

        std::cout << "\n*** Euler Implicito ***\n";
        rodaPVI2ComVariosDt("Euler Implicito",
                            [&](double dt, int nPassos, std::vector<double> &t)
                            {
                                return eulerImplicito(resolutorImplicitoPVI2, S0, T0_TAREFA, dt, nPassos, t);
                            });

        pausarPVI();
    }
}

// ----- Runge-Kutta -----
static void menuAula25()
{
    std::cout << "\n════════════════════════════════\n";
    std::cout << "   Runge-Kutta        \n";
    std::cout << "════════════════════════════════\n";
    std::cout << "  1 - RK2 ou RK3 no PVI-1 (demo)\n";
    std::cout << "  2 - PVI-2 com RK3\n";
    std::cout << "  0 - Voltar\n";

    int op = lerIntPVI("Opcao: ", 0, 2);

    if (op == 1)
    {
        int ordem = lerIntPVI("Ordem (2 ou 3): ", 2, 3);
        double dt = lerDoublePVI("Delta t: ");
        int nPassos = lerIntPVI("Numero de passos: ", 1, 100000);

        std::vector<double> tempos;
        std::vector<Estado> est;
        if (ordem == 2)
            est = rungeKutta2(F_PVI1, {2.0}, 0.0, dt, nPassos, tempos);
        else
            est = rungeKutta3(F_PVI1, {2.0}, 0.0, dt, nPassos, tempos);

        std::cout << std::fixed << std::setprecision(6);
        for (int i = 0; i < (int)tempos.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * tempos[i]);
            double err = std::fabs((est[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << tempos[i]
                      << "  y=" << est[i][0]
                      << "  y_ex=" << yEx
                      << "  erro=" << err << "%\n";
        }
        pausarPVI();
    }
    else if (op == 2)
    {
        std::cout << "\n=== Runge-Kutta 3a Ordem ===\n";
        std::cout << "PVI-2: t0=0, v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";

        Estado S0 = {V0_TAREFA, Y0_TAREFA};

        rodaPVI2ComVariosDt("RK3",
                            [&](double dt, int nPassos, std::vector<double> &t)
                            {
                                return rungeKutta3(F_PVI2, S0, T0_TAREFA, dt, nPassos, t);
                            });

        pausarPVI();
    }
}

// ----- Preditor-Corretor -----
static void menuAula26()
{
    std::cout << "\n══════════════════════════════════════════\n";
    std::cout << "    Preditor-Corretor            \n";
    std::cout << "══════════════════════════════════════════\n";
    std::cout << "  1 - Adams-Bashforth 2a ordem no PVI-1 (demo)\n";
    std::cout << "  2 - Adams-Bashforth-Moulton 3a ordem no PVI-1 (demo)\n";
    std::cout << "  3 - ABM 4a ordem no PVI-1 (demo)\n";
    std::cout << "  4 - ABM 4a ordem no PVI-2\n";
    std::cout << "  0 - Voltar\n";

    int op = lerIntPVI("Opcao: ", 0, 4);

    if (op == 1 || op == 2 || op == 3)
    {
        double dt = lerDoublePVI("Delta t: ");
        int nPassos = lerIntPVI("Numero de passos: ", 1, 100000);
        double eps = 1e-10;
        int maxIter = 50;

        std::vector<double> tempos;
        std::vector<Estado> est;

        if (op == 1)
            est = adamsBashforth2(F_PVI1, {2.0}, 0.0, dt, nPassos, eps, maxIter, tempos);
        else if (op == 2)
            est = adamsBashforthMoulton3(F_PVI1, {2.0}, 0.0, dt, nPassos, eps, maxIter, tempos);
        else
            est = adamsBashforthMoulton4(F_PVI1, {2.0}, 0.0, dt, nPassos, eps, maxIter, tempos);

        std::string nome = (op == 1) ? "AB2" : (op == 2) ? "ABM3"
                                                         : "ABM4";
        std::cout << "\n--- " << nome << " no PVI-1 ---\n";
        std::cout << std::fixed << std::setprecision(8);
        for (int i = 0; i < (int)tempos.size(); i++)
        {
            double yEx = 2.0 * std::exp((2.0 / 3.0) * tempos[i]);
            double err = std::fabs((est[i][0] - yEx) / yEx) * 100.0;
            std::cout << "  t=" << tempos[i]
                      << "  y=" << est[i][0]
                      << "  y_ex=" << yEx
                      << "  erro=" << err << "%\n";
        }
        pausarPVI();
    }
    else if (op == 4)
    {
        std::cout << "\n=== Tarefa 18 Parte 2: ABM 4a ordem no PVI-2 ===\n";
        std::cout << "PVI-2: t0=0, v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";

        Estado S0 = {V0_TAREFA, Y0_TAREFA};
        double eps = 1e-10;
        int maxIter = 50;

        rodaPVI2ComVariosDt("ABM4",
                            [&](double dt, int nPassos, std::vector<double> &t)
                            {
                                return adamsBashforthMoulton4(F_PVI2, S0, T0_TAREFA, dt, nPassos, eps, maxIter, t);
                            });

        pausarPVI();
    }
}

// Menu principal de PVI
void menuPVI()
{
    while (true)
    {
        std::cout << "\n╔══════════════════════════════════════════╗\n";
        std::cout << "║   PVI - Problemas de Valor Inicial       ║\n";
        std::cout << "╠══════════════════════════════════════════╣\n";
        std::cout << "║  1 - Euler Explicito/Implicito           ║\n";
        std::cout << "║  2 - Runge-Kutta 2a e 3a ord.            ║\n";
        std::cout << "║  3 - Preditor-Corretor                   ║\n";
        std::cout << "║  0 - Voltar                              ║\n";
        std::cout << "╚══════════════════════════════════════════╝\n";

        int op = lerIntPVI("Opcao: ", 0, 3);
        if (op == 0)
            return;
        if (op == 1)
            menuAula24();
        if (op == 2)
            menuAula25();
        if (op == 3)
            menuAula26();
    }
}