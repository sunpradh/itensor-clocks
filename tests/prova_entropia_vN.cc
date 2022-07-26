#include "itensor/all.h"
#include "clock/all.h"

using namespace itensor;

int main(int argc, char ** argv)
{
    int L = 200;
    auto sites = Clock2(L);

    Real h = 0.0;
    if (argc > 1)
        h = atof(argv[1]);

    auto H = hamiltonian(sites, {"Kinetic", -1.0, "Transv", -h, "PBC", false});

    auto sweeps = Sweeps(20);
    sweeps.maxdim() = 20, 50, 100, 200, 250, 300;
    sweeps.cutoff() = 1E-12;

    auto psi0 = randomMPS(InitState(sites));

    auto [E0, psi] = dmrg(H, psi0, sweeps, {"Silent", true});

    printfln("vN entropy for L=%d and h=%10.6f\n", L, h);
    printfln("%10s  %10s", "l", "SvN");
    println("------------------------------");
    for (auto l : range(3,L-2))
    {
        auto SvN = entropy_vN(psi, l);

        printfln("%10d  %10.7f", l, SvN);
    }
}
