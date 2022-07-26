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

int main(int argc, char ** argv)
{
    // Simulation setups
    constexpr uint N = 3; // clock order
    constexpr uint L = 20; // chain length
    auto sites = Clock<N>(L);  // clock model
    auto h_range = linspace<float>(0.0, 2.0, 20);  // coupling range
    Complex phase = exp(Complex(0, 2 * M_PI / N)); // complex phase on the long. coupling

    // Storing of the energies
    constexpr uint K = 2; // number of energy levels
    std::array<std::vector<float>, K> energies;
    for (unsigned j = 0; j < K; j++) energies[j] = std::vector<float>(h_range.size());

    // Sweeps setups
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 5, 10, 20, 50, 100;
    sweeps.cutoff() = 1E-12;

    // #pragma omp parallel for
    for (unsigned i = 0; i < h_range.size(); i++) {
        // coupling
        auto h = h_range[i];
        cout << ">>> i = " << i << ", h = " << h << "\n";

        // Hamiltonian
        auto H = hamiltonianC(sites, {
            "Kinetic",  -h,
            "Transv",   -1.0,
            "LongitRe", -h * real(1 + phase),
            "LongitIm", -h * imag(1 + phase),
            "PBC",      false
        });

        // Init state
        auto init_wf = randomMPS(InitState(sites));

        // Store the wavefunctions
        auto wfs = std::vector<MPS>();

        // Compute the ground states
        cout << "  > DMRG for ground state\n";
        auto [E0, psi0] = dmrg(H, init_wf, sweeps, {"Silent", true});
        energies[0][i] = E0;
        wfs.push_back(psi0);

        // Compute the excited states
        for (unsigned k = 1; k < K; k++) {
            cout << "  > DMRG for excited state k = " << k << "\n";
            auto [E, psi] = dmrg(H, wfs, init_wf, sweeps, {"Silent", true, "Weight", 20.0});
            energies[k][i] = E;
            wfs.push_back(psi);
        }
    }
    // Print out everything
    dump_settings::separator = ", ";
    dump(std::cout, h_range,
         energies[0],
         energies[1]
         );

    std::cout << "\nCalcoli finiti\n";

}
