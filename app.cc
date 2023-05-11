#include "itensor/all.h"
#include "clock/all.h"
#include "utils/all.h"

namespace it = itensor;
namespace ut = utils;
namespace cl = clocks;
namespace sim = clocks::simulations;

constexpr unsigned N = 3;          // clock order
// constexpr unsigned n_points = 151; // number of points to compute
constexpr unsigned n_excited = 4;  // number of excited states to compute

int main(int argc, char ** argv)
{
    // Coupling range
    constexpr unsigned n_points_full    = 151;
    constexpr unsigned n_points_focused = 201;
    auto couplings_full    = ut::linspace(1.,  1.5, n_points_full).to_array<n_points_full>();
    auto couplings_focused = ut::linspace(.75, .95, n_points_focused).to_array<n_points_focused>();

    // Chain lengths
    // auto length = 10;
    auto lengths = ut::range(50, 71);

    // Sweeps setups
    auto filename = "sweeps_input";
    auto sweeps_input = it::InputGroup(filename, "sweeps_params");
    auto sweeps_params = it::InputGroup(sweeps_input, "cutoff_based_nonoise");
    auto sweeps = it::Sweeps(7, sweeps_params);

    // Sector
    // unsigned max_sector = N/2 + 1;
    unsigned sector = 1;
    double phase_noise = 1e-4;

    for (auto length : lengths) {
        std::cout << "Clock N = " << N << ", size = " << length << "\n";
        std::cout << "------------------------------------------------------------\n";

        cout << " * Computing with no phase noise" << "\n";
        auto results_no_eps = sim::ComputeObservables<N, n_points_full, n_excited>(
                length, sweeps, couplings_full,
                {"OnlyBulk", true}
            ).compute(sector);
        results_no_eps.set_precision(12).to_csv(sim::csv_filename<N>(length, sector, "no_eps"));

        cout << " * Computing with phase noise = " << phase_noise << "\n";
        auto results_pos_eps = sim::ComputeObservables<N, n_points_focused, n_excited>(
                length, sweeps, couplings_focused,
                {"PhaseNoise", phase_noise, "OnlyBulk", true}
            ).compute(sector);
        results_pos_eps.set_precision(12).to_csv(sim::csv_filename<N>(length, sector, "eps_pos_1e4"));

        cout << " * Computing with phase noise = " << -phase_noise << "\n";
        auto results_neg_eps = sim::ComputeObservables<N, n_points_focused, n_excited>(
                length, sweeps, couplings_focused,
                {"PhaseNoise", -phase_noise, "OnlyBulk", true}
            ).compute(sector);
        results_neg_eps.set_precision(12).to_csv(sim::csv_filename<N>(length, sector, "eps_neg_1e4"));

        std::cout << "\n";
    }

    return 0;
}
