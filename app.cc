#include <string>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>

#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

using namespace itensor;
using namespace clocks;
using std::cout;

/***********************************************************/
// Simulation setups
constexpr unsigned int N = 2; // clock order
constexpr unsigned int L = 10; // chain length
constexpr unsigned int n_points = 20; // number of points to compute
constexpr unsigned int n_excited = 3; // number of excited states to compute
/***********************************************************/

using Array = std::array<double, n_points>;
using utils::Table;

template<unsigned int N>
auto dual_clock_hamiltonian(const Clock<N> & sites, double coupling, int sector) {
    Complex phase = exp(Complex(0, 2.0 * M_PI / double(N)));
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

struct Observables {
    double gs_energy;
    double disX;
    double corrZ;
};


template<unsigned int N>
Observables compute_observables(
    const Clock<N> & sites,
    const Sweeps & sweeps,
    double coupling,
    int sector
)
{
    auto H = dual_clock_hamiltonian<N>(sites, coupling, sector);
    auto init_psi = randomMPS(sites);
    auto [E0, psi] = dmrg(H, init_psi, sweeps, {"Silent", true});
    auto disX = compute_disorderC(sites, psi, "X", {1, 5}).real();
    auto corrZ = compute_correlatorC(sites, psi, "Z", "Zdag", {1, 5}).real();
    return Observables{E0, disX, corrZ};
}


template<unsigned int N>
Table<Array> compute_obs_over_range(
    const Clock<N> & sites,
    const Sweeps & sweeps,
    Array couplings,
    int sector
) {
    auto n = couplings.size();
    auto output_table = Table(
            "couplings", couplings,
            "gs_energy", Array{},
            "disX", Array{},
            "corrZ", Array{}
        );
    for (auto i : utils::range(n)) {
        auto results = compute_observables(sites, sweeps, couplings[i], sector);
        output_table["gs_energy"][i] = results.gs_energy;
        output_table["disX"][i] = results.disX;
        output_table["corrZ"][i] = results.corrZ;
    }
    return output_table;
}


string output_filename(int sector) {
    return "Z" + std::to_string(N) + "_sector_" + std::to_string(sector) + ".csv";
}



int main(int argc, char ** argv)
{
    auto sites = Clock<N>(L, {"ConserveQNs=", false});  // clock model

    // Coupling range
    Array couplings = utils::linspace(0.0, 2.0, n_points).to_array<n_points>();

    // Storage for the computed values
    Array ground_states{}, disord_X_vals{}, corr_Z_vals{};

    // Sweeps setups
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10, 20, 50, 100, 200;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;

    auto sector_results = std::array<Table<Array>, N>();

    for (auto sector : utils::range(N)) {
        cout << "Computing sector n. " << sector << "\n";
        sector_results[sector] = compute_obs_over_range(
                sites,
                sweeps,
                couplings,
                sector
            );
        cout << "\nResults:\n";
        cout << sector_results[sector].set_width(16);
        sector_results[sector].set_precision(16).to_csv(output_filename(sector));
    }

    return 0;
}
