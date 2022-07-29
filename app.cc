#include <unordered_map>
#include <vector>
#include <cmath>
#include <fstream>

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
    constexpr uint L = 3; // chain length
    auto sites = Clock<N>(L, {"ConserveQNs=", true});  // clock model

    // Hamiltonian
    float h = 2.0;
    Complex phase = exp(Complex(0, 0.4)); // complex phase on the long. coupling

    // auto H = hamiltonianC(sites, {
    //     "Kinetic",  -0,
    //     "Transv",   -1,
    //     "PBC",      false
    // });
    // PrintData(H);

    auto H = hamiltonianC(sites, {
        "KineticRe", -(1.0-h) * real(phase),
        "KineticIm", -(1.0-h) * imag(phase),
        "TransvRe",   -h * real(phase),
        "TransvIm", -h * imag(phase),
        "LongitRe", 0.0,     // -h * real(1 + phase),
        "LongitIm", 0.0,     // -h * imag(1 + phase),
        "PBC",      false
    });
    PrintData(H);

    //
    // Example of DMRG
    //

    // Sweeps setups
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 5, 10, 20, 50, 100;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;

    // Init state
    auto init_wf = randomMPS_QN(sites, 0);

    // Compute the ground states
    cout << "  > DMRG for ground state\n";
    auto [E0, psi0] = dmrg(H, init_wf, sweeps, {"Silent", true});
    cout << " E0 = " << E0 << "\n";

    return 0;
}
