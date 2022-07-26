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
    int L = 100;
    // Couplings
    Real J = 1.0; // nearest neighbour interaction
    std::vector<Real> h = range(0.0, 2.0, 0.05);

    // Set the DMRG sweeps
    auto sweeps = Sweeps(25);
    sweeps.maxdim() = 5, 10, 20, 50, 100;
    sweeps.cutoff() = 1E-8;

    // Create the chain
    auto sites = Clock<N>(L);

    // Save the GS
    auto Z = orderMPO(sites, "Z");
    auto X = orderMPO(sites, "X");

    std::vector<Real> energies;
    std::vector<Complex> Z_evs, X_evs, disZ_evs, disX_evs;

    for (auto hh : h)
    {
        // Nearest neighbour Hamiltonian
        auto H = hamiltonian(sites, {"Kinetic", -J, "Transversal", -hh, "PBC", true});
        // Initial wavefunction
        auto psi0 = randomMPS(InitState(sites));
        if (hh > 1.0)
            psi0 = InitState(sites, "0");

        // Run the DMRG algorithm
        printf("h = %6.4f ", hh);
        auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true}); print("[DMRG done] ");

        // Store the expectation values
        Z_evs.push_back(innerC(psi, Z, psi)); print("[Z done] ");
        X_evs.push_back(innerC(psi, X, psi)); print("[X done] ");
        disZ_evs.push_back(compute_disorderC(sites, psi, "Z", 5)); print("[disZ done] ");
        disX_evs.push_back(compute_disorderC(sites, psi, "X", 5)); println("[disX done] ");

        energies.push_back(energy);
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
