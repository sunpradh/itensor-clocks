#include <unordered_map>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>

// For benchmarking
#include <chrono>

#include "itensor/all.h"
#include "clock/all.h"

using namespace itensor;
using namespace clocks;
using std::cout;

int main(int argc, char ** argv)
{
    // Simulation setups
    constexpr uint N = 3; // clock order
    constexpr uint L = 10; // chain length
    auto sites = Clock<N>(L, {"ConserveQNs=", true});  // clock model

    // Hamiltonian
    float h = 2.0;
    // Complex phase = exp(Complex(0, 0.4)); // complex phase on the long. coupling

    auto H = hamiltonian(sites, {
        "Kinetic",  -1,
        "Transv",   -1.4,
        "PBC",      false
    });

    // Example of DMRG
    //----------------------------------------
    // Sweeps setups
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10, 20, 50, 100, 200;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;

    // Init state
    auto psi0 = randomMPS_QN(sites, 0);

    // Compute the ground states
    cout << "  > DMRG for ground state\n";
    auto [E0, psi] = dmrg(H, psi0, sweeps, {"Silent", true});
    cout << " E0 = " << E0 << "\n";
    //----------------------------------------

    // Example order parameters
    auto Xev = 0.5 * compute_orderC(sites, psi, "X", "Xdag");
    auto Zev = compute_orderC(sites, psi, "Z");
    cout << "<X>: " << Xev << "\n";
    cout << "<Z>: " << Zev << "\n";

    // Example disorder parameters
    auto Xdisord = compute_disorderC(sites, psi, "X", {1, 5});
    auto Zdisord = compute_disorderC(sites, psi, "Z", {1, 5});
    cout << "Xdisord: " << Xdisord << "\n";
    cout << "Zdisord: " << Zdisord << "\n";

    // Example correlators
    auto Xcorrel = compute_correlatorC(sites, psi, "X", "Xdag", {1, 5});
    auto Zcorrel = compute_correlatorC(sites, psi, "Z", "Zdag", {1, 5});
    cout << "Xcorrel: " << Xcorrel << "\n";
    cout << "Zcorrel: " << Zcorrel << "\n";

    // Example of benchmarking
    auto f = [&](){ compute_orderC(sites, psi, "X", "Xdag"); };
    auto bench = utils::benchmark<20>(f);
    bench.print_statistics<std::chrono::milliseconds>();

    return 0;
}
