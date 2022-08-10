#include <unordered_map>
#include <vector>
#include <cmath>
#include <fstream>
// For benchmarking
#include <chrono>
using namespace std::chrono;

#include "itensor/all.h"
#include "clock/all.h"
#include "clock/dump.h"

using namespace itensor;
using namespace clocks;
using std::cout;
// using std::cout;

int main(int argc, char ** argv)
{
    // Simulation setups
    constexpr uint N = 3; // clock order
    constexpr uint L = 50; // chain length
    auto sites = Clock<N>(L, {"ConserveQNs=", false});  // clock model

    // Hamiltonian
    float h = 2.0;
    // Complex phase = exp(Complex(0, 0.4)); // complex phase on the long. coupling

    auto H = hamiltonian(sites, {
        "Kinetic",  -1,
        "Transv",   -0.4,
        "PBC",      false
    });

    //
    // Example of DMRG
    //

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

    // Benchmarking
    auto f = [&](){ compute_orderC(sites, psi, "X"); };
    auto bench = benchmark<20>(f);
    bench.print_statistics<milliseconds>();

    return 0;
}
