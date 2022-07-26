#ifndef __CLOCK_CORRELATOR_H
#define __CLOCK_CORRELATOR_H

#include "clock.h"
#include "utils.h"

using Real = itensor::Real;
/************************************************************/
namespace clocks {

template<unsigned int N>
itensor::ITensor
compute_correlator_IT(
    const Clock<N> & sites,
          itensor::MPS psi,
    const string   & op_type,
    const Interval & interv
);

template<unsigned int N>
Real
compute_correlator(
    const Clock<N> & sites,
          itensor::MPS psi,
    const string   & op_type,
    const Interval & interv
);

template<unsigned int N>
Complex
compute_correlatorC(
    const clocks::Clock<N> & sites,
          itensor::MPS psi,
    const string   & op_type,
    const Interval & interv
);

}
/************************************************************/

// Compute the correlation function
// return a scalar ITensor
template<unsigned int N>
itensor::ITensor clocks::compute_correlator_IT(
    const clocks::Clock<N> & sites,
          itensor::MPS        psi,
    const string   & op_type,
    const Interval & interv
)
{
    int L = length(sites);
    auto [start, stop] = interv;

    if (start < 0 || stop > L || start >= stop)
        throw std::runtime_error("Incorrect range for disorder operator");
    if (op_type != "Z" && op_type != "X" && op_type != "Zdag" && op_type != "Xdag")
        throw std::runtime_error("Unrecognized operator for disorder operator");

    psi.position(start);

    auto Op_start = op(sites, op_type, start);
    auto Op_stop  = op(sites, op_type, stop);

    itensor::ITensor correl = psi(start);
    // Contracts the operator at the start of the interval
    correl *= Op_start;
    // find the right link index of the bra <psi|
    // Primes both the site index and the link index and then contracts
    correl *= dag(prime_inds(psi(start), "Site", rightLinkIndex(psi, start)));

    // Contracts all the site inside the interval
    for (int i=start+1; i<stop; i++)
    {
        correl *= psi(i);
        correl *= dag(prime(psi(i), "Link"));
    }

    // Contract the operators at the end of the interval
    correl *= psi(stop);
    correl *= Op_stop;
    // find the left index of the bra <psi| and then contracts with evaluated correlator
    correl *= dag(prime_inds(psi(stop), "Site", leftLinkIndex(psi, stop)));

    // Contracts all the site at the right of the interval

    return correl;
}

template<unsigned int N>
Real
clocks::compute_correlator(const clocks::Clock<N> & sites, itensor::MPS psi, const string & op_type, const Interval & interv)
{
    return elt(compute_correlator_IT(sites, psi, op_type, interv));
}

template<unsigned int N>
Complex
clocks::compute_correlatorC(const clocks::Clock<N> & sites, itensor::MPS psi, const string & op_type, const Interval & interv)
{
    return eltC(compute_correlator_IT(sites, psi, op_type, interv));
}
#endif
