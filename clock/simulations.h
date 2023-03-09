#ifndef __CLOCK_SIMULATIONS_H
#define __CLOCK_SIMULATIONS_H

#include <vector>
#include <array>
#include <string>
#include <utility>
#include <stdexcept>
#include <optional>

#include "itensor/all.h"
#include "./all.h"
#include "../utils/all.h"


namespace clocks::simulations {

namespace it = itensor;
namespace cl = clocks;
namespace ut = utils;

using complex = std::complex<double>;
using std::optional;
using std::cout;

/// Just a shorthand
template<unsigned N>
    using Clock = cl::Clock<N>;

/// Dynamic size array

/// Quick conversion to string
template<typename T>
string str(T a) {
    return std::to_string(a);
}

/// Filename for the output csv
template<unsigned int N>
string output_csv_filename(unsigned len, unsigned sector) {
    return "Z" + str(N) + "_L_" + str(len) + "_sector_" + str(sector) + ".csv";
}

template<unsigned N, unsigned n_points, unsigned n_excited = 1>
struct Compute {

    using Array = std::array<double, n_points>;
    using Vector = std::vector<double>;
    using Table = ut::Table<Array>;

    unsigned size;
    Clock<N> sites;
    Array couplings;
    it::Sweeps sweeps;
    it::Args args;

    /// Constructor
    /// needs chain length, sweeps and couplings
    Compute(
        unsigned chain_length_,
        const it::Sweeps & sweeps_,
        const Array & couplings_,
        const it::Args & args_ = {}
    ) :
        size(chain_length_),
        sites(chain_length_, {"ConserveQNs", false}),
        couplings(couplings_),
        sweeps(sweeps_),
        args(args_) {};

    /// Struct for storing the result of a single DMRG calculation
    struct Observables {
        double gs_energy;
        optional<double> disorder;
        optional<double> correlator_half;
        optional<Vector> correlator;
        optional<Vector> excited_energies;
    };

    /// Dual Clock Hamiltonian
    auto dual_hamiltonian(double coupling, unsigned sector) {
        complex phase = exp(complex(0.0, 2.0 * M_PI * double(sector) / double(N)));
        complex longit_factor = 1.0 + phase;
        return hamiltonianC<N>(sites, {
                "Kinetic",  - coupling,
                "Transv",   - 1.0,
                "LongitRe", - coupling * longit_factor.real(),
                "LongitIm", - coupling * longit_factor.imag(),
                "PBC",      args.getBool("PBC", false)
            });
    }

    /// Disorder operator, equivalent to the Wilson loop
    optional<double> disorder(it::MPS & psi) {
        if (args.getBool("NoDisorder", false))
            return {};
        return 0.5 * abs(cl::compute_disorderC(sites, psi, "X", {1, size/2}) +
                         cl::compute_disorderC(sites, psi, "Xdag", {1, size/2}));
    };

    /// Correlator from the start to the middle of the chain, equivalent to the 't Hooft string
    optional<double> half_chain_correlator(it::MPS & psi) {
        if (args.getBool("NoHalfChainCorrelator", false))
            return {};
        return abs(cl::compute_correlatorC(sites, psi, "Z", "Zdag", {1, size/2}));
    };

    /// Compute correlator on a given range inside the chain
    optional<Vector> correlator(it::MPS & psi, unsigned begin, unsigned end){
        if (args.getBool("NoCorrelator", false))
            return {};
        if (begin >= end)
            throw std::invalid_argument("`begin` must be strictly smaller than `end`");

        Vector corr_values{};
        corr_values.reserve(end - begin);
        for (auto pos : ut::range(begin+1, end))
            corr_values.emplace_back(
                    abs(cl::compute_correlatorC(sites, psi, "Z", "Zdag", {begin, pos}))
                );
        return corr_values;
    };

    /// Compute energy of the excited levels
    optional<Vector> excited_levels(it::MPO & hamiltonian, it::MPS & psi0, it::MPS & init_psi) {
        if (args.getBool("NoExcitedLevels", false))
            return {};
        Vector excited_energies(n_excited);
        auto wavefunctions = std::vector<it::MPS>{};
        wavefunctions.reserve(n_excited);
        wavefunctions.push_back(psi0);

        // for (auto& level : excited_energies) {
        for (unsigned n=0; n < n_excited; n++) {
            auto [E, psi] = dmrg(
                    hamiltonian,
                    wavefunctions,
                    init_psi,
                    sweeps,
                    {"Silent", true, "Weight", 10.0}
                );
            excited_energies.at(n) = E;
            wavefunctions.push_back(psi);
        }

        return excited_energies;
    };

    /// Compute the observables for a given coupling and sector
    std::pair<Observables, it::MPS> observables_at(
        it::MPS & init_psi,
        double coupling,
        unsigned sector
    ) {
        auto H = dual_hamiltonian(coupling, sector);
        auto [gs_energy, psi] = dmrg(H, init_psi, sweeps, {"Silent", true});
        auto results = Observables{
            gs_energy,
            disorder(psi),
            half_chain_correlator(psi),
            correlator(psi, size/4, 3*size/4),
            excited_levels(H, psi, init_psi)
        };

        return std::make_pair(results, psi);
    };

    /// Computing observables for each couplings for a given sector
    Table all(unsigned sector) {
        unsigned corr_begin = size/4, corr_end = 3*size/4;
        auto results = new_table(corr_begin, corr_end);
        auto init_psi = it::randomMPS(sites);

        std::cout << " * Computing sector n. " << sector << "\n";
        unsigned n_steps = couplings.size(), step = 0;
        // DMRG calculation for each coupling
        #pragma omp parallel for
        for (auto i : ut::range(n_steps)) {
            #pragma omp critical
            { print_progress(++step, n_steps); }
            auto [obs, psi] = observables_at(init_psi, couplings.at(i), sector);
            fill_table_row(results, obs, i);
            init_psi = psi;
        }
        std::cout << " Done!\n";

        return results;
    };

private:
    /// Create the Table object for storing the results of a single
    /// DMRG calculation
    Table new_table(unsigned corr_begin, unsigned corr_end) {
        auto table = Table("couplings", couplings, "gs_energy", Array{});

        // Optional columns
        if (!args.getBool("NoDisorder", false))
            table.add_columns("disorder", Array{});
        if (!args.getBool("NoHalfChainCorrelator", false))
            table.add_columns("corr_half", Array{});
        // Optional excited energies columns
        if (!args.getBool("NoExcitedLevels", false))
            for (auto n : ut::range(n_excited))
                table.add_columns("E" + str(n+1), Array{});
        // Optional correlator columns
        if (!args.getBool("NoCorrelator", false))
            for (auto r : ut::range(1u, corr_end - corr_begin))
                table.add_columns("corr_R_" + str(r), Array{});
        return table;
    }

    /// Fill all the correlator entries of the given row
    void fill_correlator_row(Table & table, const Vector & corr_values, unsigned row) {
        for (auto r : ut::range(corr_values.size()))
           table["corr_R_" + str(r+1)][row] = corr_values.at(r);
    }

    /// Fills all the excited energies entries of a given row
    void fill_excited_row(Table & table, const Vector & excited_levels, unsigned row) {
        for (auto n : ut::range(n_excited))
            table["E" + str(n+1)][row] = excited_levels.at(n);
    }

    /// Fill all the entries of a row with the given observables
    void fill_table_row(Table & table, const Observables & obs, unsigned row) {
        table["gs_energy"][row] = obs.gs_energy;
        if (obs.disorder)
            table["disorder"][row] = *obs.disorder;
        if (obs.correlator_half)
            table["corr_half"][row] = *obs.correlator_half;
        if (obs.excited_energies)
            fill_excited_row(table, *obs.excited_energies, row);
        if (obs.correlator)
            fill_correlator_row(table, *obs.correlator, row);
    }

    /// Simply print the progress
    void print_progress(unsigned n, unsigned total) {
        std::cout << "\033[2K\r"
            << "   In progress [" << n << "/" << total << "]"
            << std::flush;
    }
};

}

#endif
