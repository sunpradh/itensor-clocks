#include <string>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>

#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

// Simulation setups
//----------------------------------------------------------

constexpr unsigned int N = 3; // clock order
constexpr unsigned int n_points = 51; // number of points to compute
constexpr unsigned int n_excited = 2; // number of excited states to compute

// Types and namespaces
//----------------------------------------------------------

namespace it = itensor;
namespace cl = clocks;
namespace ut = utils;

using std::complex;
using std::cout;
using ut::Table;
template<unsigned int N> using Clock = cl::Clock<N>;
using Array = std::array<double, n_points>;
using EnergyArray = std::array<double, n_excited>;


struct Observables {
    double gs_energy;
    double disX;
    double corrZ;
    EnergyArray excited_energies;
};


// Functions
//----------------------------------------------------------


// Quick conversion to string
template<typename T>
string str(T a) {
    return std::to_string(a);
}

// Filename for the output csv
string output_filename(int len, int sector) {
    return "Z" + str(N) + "_L_" + str(len) + "_sector_" + str(sector) + ".csv";
}

// Clock Hamiltonian dual to the ladder LGT
template<unsigned int N>
auto dual_clock_hamiltonian(const Clock<N> & sites, double coupling, int sector) {
    complex phase = exp(complex(0.0, 2.0 * M_PI * double(sector) / double(N)));
    complex longit_factor = 1.0 + phase;

    return hamiltonianC<N>(
        sites,
        {
            "Kinetic",  - coupling,
            "Transv",   - 1.0,
            "LongitRe", - coupling * longit_factor.real(),
            "LongitIm", - coupling * longit_factor.imag()
        }
    );
}

// Compute excited energies (above the ground state)
template<unsigned int N>
EnergyArray compute_excited_energies(
        const Clock<N> & sites,
        it::MPO & hamiltonian,
        it::MPS & psi0,
        const it::Sweeps & sweeps
) {
    EnergyArray excited_energies{0};
    auto wavefunctions = std::vector<MPS>{psi0};

    for (auto& level : excited_energies) {
        auto [E, psi] = dmrg(
                hamiltonian,
                wavefunctions,
                randomMPS(sites),
                sweeps,
                {"Silent", true}
            );
        level = E;
        wavefunctions.push_back(psi);
    }

    return excited_energies;
}

// Compute all the required observables (see members of Observables struct)
template<unsigned int N>
Observables compute_observables(
    const Clock<N> & sites,
    const it::Sweeps & sweeps,
    double coupling,
    int sector
) {
    auto len = length(sites);
    auto H = dual_clock_hamiltonian<N>(sites, coupling, sector);
    auto init_psi = randomMPS(sites);
    auto [E0, psi] = dmrg(H, init_psi, sweeps, {"Silent", true});
    auto disX = 0.5 * abs(
            compute_disorderC(sites, psi, "X", {1, int(len/2)}) +
            compute_disorderC(sites, psi, "Xdag", {1, int(len/2)})
        );
    auto corrZ = abs(compute_correlatorC(sites, psi, "Z", "Zdag", {1, int(len/2)}));
    auto excited_energies = compute_excited_energies(sites, H, psi, sweeps);
    return Observables{E0, disX, corrZ, excited_energies};
}


// Compute the observables over a given range of couplings
template<unsigned int N>
Table<Array> compute_obs_over_range(
    const Clock<N> & sites,
    const it::Sweeps & sweeps,
    const Array & couplings,
    int sector
) {
    auto output = Table(
            "couplings", couplings,
            "gs_energy", Array{},
            "disX", Array{},
            "corrZ", Array{}
        );
    for (auto n : utils::range(n_excited))
        output.add_columns("E" + str(n+1), Array{});

    auto n_steps = couplings.size();
    int n = 0;
    #pragma omp parallel for
    for (auto i : utils::range(n_steps)) {
        #pragma omp critical
        {
            n++;
            cout << "\033[2K\r"
                << "   In progress [" << n << "/" << n_steps << "]"
                << std::flush;
        }

        auto results = compute_observables(sites, sweeps, couplings.at(i), sector);

        output["gs_energy"][i] = results.gs_energy;
        output["disX"][i]      = results.disX;
        output["corrZ"][i]     = results.corrZ;
        for (auto n : utils::range(n_excited))
            output["E" + str(n+1)][i] = results.excited_energies.at(n);
    }
    cout << " Done!\n";

    return output;
}

template<unsigned int N>
auto compute_obs(
    const Clock<N> & sites,
    const it::Sweeps & sweeps,
    const Array & couplings
) {
    constexpr unsigned int max_sector = N/2 + 1;
    std::array<Table<Array>, max_sector> results{};

    for (auto sector : utils::range(max_sector)) {
        cout << " * Computing sector n. " << sector << "\n";
        results[sector] = compute_obs_over_range(
                sites,
                sweeps,
                couplings,
                sector
            );

        // cout << "\nResults:\n" << results[sector].set_width(16) << "\n";
        results[sector].set_precision(16).to_csv(output_filename(length(sites), sector));
    }

    return results;
}


// Main function
//----------------------------------------------------------

int main(int argc, char ** argv)
{
    // Sweeps setups
    auto sweeps = it::Sweeps(5);
    sweeps.maxdim() = 10, 20, 50, 100, 150;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;

    // Coupling range
    Array couplings = utils::linspace(0.0, 2.0, n_points).to_array<n_points>();

    // Chain lengths
    auto lengths = {10, 20, 30, 40};
    // auto lengths = {10};

    for (auto L : lengths) {
        // Clock model (no QN conservation)
        cout << "\nClock model N=" << N << ", chain size = " << L << "\n";
        cout << "--------------------------------------------------\n";
        auto sites = Clock<N>(L, {"ConserveQNs=", false});
        compute_obs(sites, sweeps, couplings);
    }

    return 0;
}
