#ifndef __CLOCK_CORRELATOR_H
#define __CLOCK_CORRELATOR_H

#include "clock.h"
#include "../utils/all.h"

using Real = itensor::Real;
using itensor::MPS;
using itensor::ITensor;
/************************************************************/
namespace clocks {

template<unsigned int N>
ITensor
compute_correlator_IT(
    const Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
);

template<unsigned int N>
Real
compute_correlator(
    const Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
);

template<unsigned int N>
Complex
compute_correlatorC(
    const Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
);

}
/************************************************************/

// Compute the correlation function
// return a scalar ITensor
template<unsigned int N>
ITensor clocks::compute_correlator_IT(
    const clocks::Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
)
{
    if (!is_valid_op(op1))
        throw std::runtime_error("Unrecognized first operator for correlator");
    if (!is_valid_op(op2))
        throw std::runtime_error("Unrecognized second operator for correlator");

    int L = length(sites);
    auto [begin, end] = interv;
    if (begin < 0 || end > L || begin >= end)
        throw std::runtime_error("Incorrect position for correlator");

    auto op_first  = op(sites, op1, begin);
    auto op_second = op(sites, op2, end);

    psi.position(begin);
    ITensor correl = psi(begin);
    // Contracts the operator at the start of the interval
    correl *= op_first;

    // find the right link index of the bra <psi|
    // Primes both the site index and the link index and then contracts
    correl *= dag(utils::prime_inds(psi(begin), "Site", rightLinkIndex(psi, begin)));

    // Contracts all the site inside the interval
    for (int i=begin+1; i<end; i++)
        correl *= psi(i) * dag(prime(psi(i), "Link"));

    // Contract the operators at the end of the interval
    correl *= psi(end);
    correl *= op_second;

    // find the left index of the bra <psi| and then contracts with evaluated correlator
    correl *= dag(utils::prime_inds(psi(end), "Site", leftLinkIndex(psi, end)));
    // Contracts all the site at the right of the interval

    return correl;
}

template<unsigned int N>
Real
clocks::compute_correlator(
    const Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
)
{
    return elt(compute_correlator_IT(sites, psi, op1, op2, interv));
}

template<unsigned int N>
Complex
clocks::compute_correlatorC(
    const Clock<N> & sites,
          MPS      & psi,
    const string   & op1,
    const string   & op2,
    const Interval & interv
)
{
    return eltC(compute_correlator_IT(sites, psi, op1, op2, interv));
}
#endif
