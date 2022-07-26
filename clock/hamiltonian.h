#ifndef __CLOCK_HAMILTONIAN_H
#define __CLOCK_HAMILTONIAN_H

#include "itensor/all.h"
#include "clock.h"

using Real = itensor::Real;
/************************************************************/
namespace clocks {

// Non-chiral Hamiltonians (aka only real couplings).
// pass parameters explicitly
template<unsigned int N>
itensor::MPO
hamiltonian(
    const Clock<N> & sites,
    Real kinet  = 1.0,
    Real transv = 0.0,
    Real longit = 0.0,
    bool pbc    = false
);

// Non-chiral Hamiltonians (aka only real couplings).
// pass Args object for parameters
template<unsigned int N>
itensor::MPO
hamiltonian(
    const clocks::Clock<N> & sites,
    const itensor::Args & args = itensor::Args::global()
);

// Chiral Hamiltonians (admits complex couplings)
// Only for N>=3
template<unsigned int N>
itensor::MPO
hamiltonianC(
    const clocks::Clock<N> & sites,
    Complex kinet  = Complex(-1, 0),
    Complex transv = Complex( 0, 0),
    Complex longit = Complex( 0, 0),
    bool    pbc    = false
);

template<unsigned int N>
itensor::MPO
hamiltonianC(
    const clocks::Clock<N> & sites,
    const itensor::Args & args = itensor::Args::global()
);

}
/************************************************************/


template<unsigned int N>
itensor::MPO
clocks::hamiltonian(const clocks::Clock<N> & sites, Real kin, Real transv, Real longit, bool pbc)
{
    int L = length(sites);
    auto H_ampo = itensor::AutoMPO(sites);

    // Kinetic term
    for (int i=1; i < L; i++)
    {
        H_ampo += kin, "Zdag", i+1, "Z",    i;
        H_ampo += kin, "Z",    i+1, "Zdag", i;
    }
    // naive approach to PBC for the kinetic term
    if (pbc)
    {
        H_ampo += kin, "Zdag", L, "Z", 1;
        H_ampo += kin, "Z", L, "Zdag", 1;
    }

    // Transversal field
    if (transv != 0.0)
        for (int i=1; i <= L; i++)
        {
            H_ampo += transv, "X",    i;
            H_ampo += transv, "Xdag", i;
        }

    // Longitudinal field
    if (longit != 0.0)
        for (int i=1; i<=L; i++)
        {
            H_ampo += longit, "Z", i;
            H_ampo += longit, "Zdag", i;
        }

    return toMPO(H_ampo);
}

template<unsigned int N>
itensor::MPO
clocks::hamiltonian(const clocks::Clock<N> & sites, const itensor::Args & args)
{
    auto pbc    = args.getBool("PBC",     false);
    auto kin    = args.getReal("Kinetic", -1.0);
    auto transv = args.getReal("Transv",   0.0);
    auto longit = args.getReal("Longit",   0.0);

    return hamiltonian(sites, kin, transv, longit, pbc);
}

template<unsigned int N>
itensor::MPO
clocks::hamiltonianC(const clocks::Clock<N> & sites, Complex kin, Complex transv, Complex longit, bool pbc)
{
    int L = itensor::length(sites);
    auto H_ampo = itensor::AutoMPO(sites);

    // Kinetic term
    for (int i=1; i < L; i++)
    {
        H_ampo += kin,       "Zdag", i+1, "Z",    i;
        H_ampo += conj(kin), "Z",    i+1, "Zdag", i;
    }
    // naive approach to PBC for the kinetic term
    if (pbc)
    {
        H_ampo += kin,       "Zdag", L, "Z",    1;
        H_ampo += conj(kin), "Z",    L, "Zdag", 1;
    }

    // Transversal field
    if (transv != Complex(0.0))
        for (int i=1; i<=L; i++)
        {
            H_ampo += transv,       "X",    i;
            H_ampo += conj(transv), "Xdag", i;
        }

    // Longitudinal field
    if (longit != Complex(0.0))
        for (int i=1; i<=L; i++)
        {
            H_ampo += longit,       "Z",    i;
            H_ampo += conj(longit), "Zdag", i;
        }

    return toMPO(H_ampo);
}

// Chiral Hamiltonian (only for N>=3)
template<unsigned int N>
itensor::MPO
clocks::hamiltonianC(const clocks::Clock<N> & sites, const itensor::Args& args)
{
    auto pbc      = args.getBool("PBC",       false);
    auto kinRe    = args.getReal("KineticRe", args.getReal("Kinetic", -1.0));
    auto transvRe = args.getReal("TransvRe",  args.getReal("Transv",  -1.0));
    auto longitRe = args.getReal("LongitRe",  args.getReal("Longit",  -1.0));
    auto kinIm    = args.getReal("KineticIm", 0.0);
    auto transvIm = args.getReal("TransvIm",  0.0);
    auto longitIm = args.getReal("LongitIm",  0.0);

    return hamiltonianC(
                sites,
                Complex(kinRe,    kinIm),
                Complex(transvRe, transvIm),
                Complex(longitRe, longitIm),
                pbc
            );
}

#endif
