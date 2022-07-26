#include "itensor/all.h"
#include "clock/all.h"
#include <unordered_map>
#include <cmath>

template<typename T1, typename T2> using map = std::unordered_map<T1, T2>;

/*
* Valuta quanti sweeps sono necessari affinchè l'entropia di entanglement non varii
* di un certo cutoff, al variare della taglia della catena
*
* Può essere riciclato in realtà per qualsiasi osservabile
*/
void optimal_sweeps()
{
    Real EE_cutoff = 1E-10;

    auto Lrange = range(40,201,20);
    auto sweeprange = range(5,41);
    map<int,Real> sweeps_needed;

    print("progress: ");
    #pragma omp parallel for
    for (auto L : Lrange)
    {
        auto sites = Clock<2>(L);
        auto H = hamiltonian(sites, {"Transv", -1.0});
        auto psi0 = randomMPS(InitState(sites));
        Real EEold = -1, EE = 0;
        int nsweep = 0;
        for (auto s : sweeprange)
        {
            auto sweeps = Sweeps(s);
            sweeps.maxdim() = 10, 20, 50, 100, 200, 300;
            sweeps.cutoff() = 1E-10;
            auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true}) ;
            EE = entropy_vN(psi, int(L/2));
            if ( abs(EE - EEold) > EE_cutoff )
                EEold = EE;
            else
            {
                nsweep = s;
                break;
            }
        }
        if (nsweep == -1) nsweep = sweeprange.back();
        sweeps_needed.insert({L, nsweep});
        print("#");
    }
    println(" done\n");

    printfln("%5s\t%5s", "L", "sweep");
    for (auto L : Lrange)
        printfln("%5d\t%5d", L, sweeps_needed.at(L));
}

/*
* Misura l'entropia di entanglement al variare del numero di sweep
* fissata una certa taglia della catena e al punto critico h=1
*/
void EE_vs_sweep()
{
    int L=200;
    auto sites = Clock<2>(L);
    auto H = hamiltonian(sites, {"Transv", -1});
    auto psi0 = randomMPS(InitState(sites));
    std::unordered_map<int,Real> entropies;

    auto nsweeps = range(5,41,5);
    print("progress: ");
    #pragma omp parallel for
    for (auto s : nsweeps)
    {
        auto sweeps = Sweeps(s);
        sweeps.maxdim() = 10, 20, 50, 100, 200;
        sweeps.cutoff() = 1E-10;
        auto [E0, psi] = dmrg(H, psi0, sweeps, {"Silent", true});
        entropies.insert({s, entropy_vN(psi, int(L/2))});
        print('#');
    }
    println(" done");
    for (auto s : nsweeps)
        printfln("%5d\t%15.12f", s, entropies.at(s));
}


int main()
{
    optimal_sweeps();
}

