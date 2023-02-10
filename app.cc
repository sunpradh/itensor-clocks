#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <utility>


#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

using namespace itensor;
using namespace clocks;
using std::cout;


template<int N>
auto dual_clock_hamiltonian(const Clock<N> & sites, float coupling, int sector) {
    Complex phase = exp(Complex(0, 2.0 * M_PI / float(N)));
    Complex longit_factor = (1.0 + pow(phase, sector));

    return hamiltonianC(
        sites,
        {
            "Kinetic", - coupling,
            "Transv", -1.0,
            "LongitRe", - coupling * longit_factor.real(),
            "LongitIm", - longit_factor.imag()
        }
    );
}


int main(int argc, char ** argv)
{
    // Simulation setups
    constexpr uint N = 2; // clock order
    constexpr uint L = 10; // chain length
    constexpr uint n_points = 20; // number of points to compute
    constexpr uint n_excited = 3; // number of excited states to compute

    using Array = std::array<float, n_points>;
    auto sites = Clock<N>(L, {"ConserveQNs=", false});  // clock model

    // Coupling range
    auto coupling_range = utils::linspace(0.0f, 2.0f, n_points);

    auto couplings = Array{};
    for (auto i : utils::range<int>(n_points)) couplings[i] = coupling_range[i];

    // Storage for the computed values
    Array ground_states{}, disord_X_vals{}, corr_Z_vals{};

    // Sweeps setups
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10, 20, 50, 100, 200;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;

    for (auto i : utils::range<int>(n_points)) {
        // Hamiltonian
        auto H = dual_clock_hamiltonian<N>(sites, couplings[i], 0);

        // Random init state
        auto psi0 = randomMPS(sites);
        // Compute the ground states
        auto [E0, psi] = dmrg(H, psi0, sweeps, {"Silent", true});
        ground_states[i] = E0;

        // disorder operator
        disord_X_vals[i] = compute_disorderC(sites, psi, "X", {1, 5}).real();
        // correlators
        corr_Z_vals[i] = compute_correlatorC(sites, psi, "Z", "Zdag", {1, 5}).real();
    }

    utils::Table table(
            "couplings", couplings,
            "E0",        ground_states,
            "disX",      disord_X_vals,
            "corrZ",     corr_Z_vals
        );
    table.set_width(16).print();
    table.set_precision(16).to_csv("out.csv");

    return 0;
}
