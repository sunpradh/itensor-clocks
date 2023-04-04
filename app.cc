#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

namespace it = itensor;
namespace ut = utils;
namespace cl = clocks;
namespace sim = clocks::simulations;

constexpr unsigned N = 3;          // clock order
constexpr unsigned n_points = 151; // number of points to compute
constexpr unsigned n_excited = 1;  // number of excited states to compute

int main(int argc, char ** argv)
{
    // Coupling range
    auto couplings = ut::linspace(.0, 1., n_points).to_array<n_points>();

    // Chain lengths
    // auto length = 10;
    auto lengths = ut::range(12, 31, 2);

    // Sweeps setups
    auto filename = "sweeps_input";
    auto sweeps_input = it::InputGroup(filename, "sweeps_params");
    auto sweeps_params = it::InputGroup(sweeps_input, "cutoff_based_nonoise");
    auto sweeps = it::Sweeps(7, sweeps_params);

    // Sector
    // unsigned max_sector = N/2 + 1;
    unsigned sector = 1;

    for (auto length : lengths) {
        std::cout << "Clock N=" << N << ", size = " << length << "\n";
        std::cout << "------------------------------------------------------------\n";
        auto compute = sim::ComputeObservables<N, n_points, n_excited>(
                length,
                sweeps,
                couplings,
                {"NoExcited", true}
            );
        auto results = compute.all(sector);
        results.set_precision(12).to_csv(sim::csv_filename<N>(length, sector));
        std::cout << "\n";
    }

    return 0;
}
