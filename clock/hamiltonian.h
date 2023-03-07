#ifndef __CLOCK_HAMILTONIAN_H
#define __CLOCK_HAMILTONIAN_H

#include "itensor/all.h"
#include "clock.h"

/************************************************************/
namespace clocks {

// Non-chiral Hamiltonians (aka only real couplings).
// pass parameters explicitly
template<unsigned int N>
it::MPO hamiltonian(
    const Clock<N> & sites,
    double kinet  = 1.0,
    double transv = 0.0,
    double longit = 0.0,
    bool pbc    = false
);

/// Non-chiral Hamiltonians (aka only real couplings).
/// pass Args object for parameters
template<unsigned int N>
it::MPO hamiltonian(
    const clocks::Clock<N> & sites,
    const it::Args & args = it::Args::global()
);

/// Chiral Hamiltonians (admits complex couplings)
/// Only for N>=3, pass parameters explicitly
template<unsigned int N>
it::MPO hamiltonianC(
    const clocks::Clock<N> & sites,
    complex kinet  = complex(-1, 0),
    complex transv = complex( 0, 0),
    complex longit = complex( 0, 0),
    bool    pbc    = false
);

/// Chiral Hamiltonians (admits complex couplings)
/// Only for N>=3, pass Args for parameters
template<unsigned int N>
it::MPO hamiltonianC(
    const clocks::Clock<N> & sites,
    const it::Args & args = it::Args::global()
);

/************************************************************/


//
// Real Hamiltonian
//
template<unsigned int N>
it::MPO hamiltonian(
    const clocks::Clock<N> & sites,
    double kin,
    double transv,
    double longit,
    bool pbc
) {
    int L = length(sites);
    auto H_ampo = it::AutoMPO(sites);

    // Kinetic term
    for (int i=1; i < L; i++) {
        H_ampo += kin, "Zdag", i+1, "Z",    i;
        H_ampo += kin, "Zdag", i,   "Z", i+1;
    }
    // naive approach to PBC for the kinetic term
    if (pbc) {
        H_ampo += kin, "Zdag", L, "Z", 1;
        H_ampo += kin, "Zdag", 1, "Z", L;
    }

    // Transversal field
    if (transv != 0.0)
        for (int i=1; i <= L; i++) {
            H_ampo += transv, "X",    i;
            H_ampo += transv, "Xdag", i;
        }

    // Longitudinal field
    if (longit != 0.0)
        for (int i=1; i<=L; i++) {
            H_ampo += longit, "Z", i;
            H_ampo += longit, "Zdag", i;
        }

    return toMPO(H_ampo);
}


template<unsigned int N>
it::MPO hamiltonian(
    const clocks::Clock<N> & sites,
    const it::Args    & args
) {
    auto pbc    = args.getBool("PBC",     false);
    auto kin    = args.getReal("Kinetic", -1.0);
    auto transv = args.getReal("Transv",   0.0);
    auto longit = args.getReal("Longit",   0.0);

    return hamiltonian(sites, kin, transv, longit, pbc);
}

//
// Chiral Hamiltonian (only for N>=3),
// i.e. Hamiltonian with complex coeffiecents
//
template<unsigned int N>
it::MPO hamiltonianC(
    const clocks::Clock<N> & sites,
    complex kin,
    complex transv,
    complex longit,
    bool pbc
) {
    int L = it::length(sites);
    auto H_ampo = it::AutoMPO(sites);

    // Kinetic term
    for (int i=1; i < L; i++) {
        H_ampo += kin,       "Zdag", i+1, "Z", i;
        H_ampo += conj(kin), "Zdag", i,   "Z", i+1;
    }
    // naive approach to PBC for the kinetic term
    if (pbc) {
        H_ampo += kin,       "Zdag", L, "Z", 1;
        H_ampo += conj(kin), "Zdag", 1, "Z", L;
    }

    // Transversal field
    if (transv != complex(0.0))
        for (int i=1; i<=L; i++) {
            H_ampo += transv,       "X",    i;
            H_ampo += conj(transv), "Xdag", i;
        }

    // Longitudinal field
    if (longit != complex(0.0))
        for (int i=1; i<=L; i++) {
            H_ampo += longit,       "Z",    i;
            H_ampo += conj(longit), "Zdag", i;
        }

    return toMPO(H_ampo);
}


template<unsigned int N>
it::MPO hamiltonianC(
    const clocks::Clock<N> & sites,
    const it::Args    & args
) {
    auto pbc      = args.getBool("PBC",       false);
    auto kinRe    = args.getReal("KineticRe", args.getReal("Kinetic", -1.0));
    auto transvRe = args.getReal("TransvRe",  args.getReal("Transv",  -1.0));
    auto longitRe = args.getReal("LongitRe",  args.getReal("Longit",  -1.0));
    auto kinIm    = args.getReal("KineticIm", 0.0);
    auto transvIm = args.getReal("TransvIm",  0.0);
    auto longitIm = args.getReal("LongitIm",  0.0);

    return hamiltonianC(
        sites,
        complex(kinRe,    kinIm),
        complex(transvRe, transvIm),
        complex(longitRe, longitIm),
        pbc
    );
}

}

#endif
