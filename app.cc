#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

namespace it = itensor;
namespace ut = utils;
namespace cl = clocks;
namespace sim = clocks::simulations;

constexpr unsigned N = 2; // clock order
constexpr unsigned n_points = 50; // number of points to compute
constexpr unsigned n_excited = 2; // number of excited states to compute

int main(int argc, char ** argv)
{
    // Sweeps setups
    auto sweeps = it::Sweeps(5);
    sweeps.maxdim() = 10, 20, 50, 100, 150;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 4;
    sweeps.noise() = 1e-7, 1e-8, 0.0;

    // Coupling range
    auto couplings = ut::linspace(0.0, 2.0, n_points).to_array<n_points>();

    // Chain lengths
    auto lengths = {10, 20, 30};

    unsigned max_sector = N/2 + 1;
    for (auto L : lengths) {
        std::cout << "Clock, size = " << L << "\n";
        std::cout << "------------------------------------------------------------\n";
        for (auto sector : ut::range(max_sector)) {
            auto results = sim::Compute<N, n_points, n_excited>(L, sweeps, couplings).all(sector);
            results.set_precision(12).to_csv(sim::output_csv_filename<N>(L, sector));
        }
        std::cout << "\n";
    }

    return 0;
}
