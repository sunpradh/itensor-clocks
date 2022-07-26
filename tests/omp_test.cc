/*
 * Quick test per vedere se clock.h funziona oppure no
 */
#include <vector>
#include <string>
#include <stdexcept>
#include "../clock/all.h"
#include "../clock/utils.h"
#include "utils.h"

int main()
{
    // Order of the clock
    constexpr unsigned int N = 2;
    // Lenght of the chain
    int L = 200;
    // Couplings
    Real J = 1.0; // nearest neighbour interaction
    std::vector<Real> h = range(0.0, 2.0, 0.05);

    // Set the DMRG sweeps
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 5, 10, 20, 50, 100;
    sweeps.cutoff() = 1E-8;

    // Create the chain
    auto sites = Clock<N>(L);

    // Save the GS
    auto Z = orderMPO(sites, "Z");
    auto X = orderMPO(sites, "X");

    std::vector<Real> energies(h.size());
    std::vector<Complex> Z_evs(h.size()),
                         X_evs(h.size()),
                         disZ_evs(h.size()),
                         disX_evs(h.size());

    // Parallelize with OMP
    // the number of threads can be governed with the env variable OMP_NUM_THREADS
    #pragma omp parallel for
    for (unsigned int i=0; i < h.size(); i++)
    {
        // Nearest neighbour Hamiltonian
        auto H = hamiltonian(sites, {"Kinetic", -J, "Transversal", -h[i], "PBC", false});
        // Initial wavefunction
        auto psi0 = randomMPS(InitState(sites));
        if (h[i] > 1.0)
            psi0 = InitState(sites, "0");

        // Run the DMRG algorithm
        printfln("h = %6.4f ", h[i]);
        auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true});

        // Store the expectation values
        energies[i] = energy;
        Z_evs[i]    = innerC(psi, Z, psi);
        X_evs[i]    = innerC(psi, X, psi);
        disZ_evs[i] = compute_disorderC(sites, psi, "Z", 5);
        disX_evs[i] = compute_disorderC(sites, psi, "X", 5);
    }

    println("\n RESULTS \n");
    printfln("%10s %12s %10s %10s %10s %10s", "h", "E0", "Z", "X", "disZ", "disX");
    int h_size = h.size();
    for (int n=0; n < h_size; n++)
        printfln("%10.7f %12.6f %10.7f %10.7f %10.7f %10.7f",
                h[n],
                energies[n],
                abs(Z_evs[n]),
                abs(X_evs[n]),
                abs(disZ_evs[n]),
                abs(disX_evs[n])
            );

}
