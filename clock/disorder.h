#ifndef __CLOCK_DISORDER_H
#define __CLOCK_DISORDER_H

#include "clock.h"
#include "utils.h"

/************************************************************/

namespace clocks {

template<unsigned int N>
itensor::ITensor
compute_disorder_IT(
    const Clock<N> & sites,
    itensor::MPS        psi,
    const string   & op_type,
    const Interval & interv
);


template<unsigned int N>
Real
compute_disorder(
    const Clock<N> & sites,
    itensor::MPS        psi,
    const string   & op_type,
    const Interval & interv
);

template<unsigned int N>
Complex
compute_disorderC(
    const Clock<N> & sites,
    itensor::MPS        psi,
    const string   & op_type,
    const Interval & interv
);

}
/************************************************************/

// Compute disorder parameter
// Return a scalar ITensor
template<unsigned int N>
itensor::ITensor
clocks::compute_disorder_IT(const clocks::Clock<N> & sites, itensor::MPS psi, const string & op_type, const Interval & interv)
{
    int L = length(sites);
    auto [start, stop] = interv;

    if (start < 0 || stop > L || start >= stop)
        throw std::runtime_error("Incorrect range for disorder operator");
    if (op_type != "Z" && op_type != "X" && op_type != "Zdag" && op_type != "Xdag")
        throw std::runtime_error("Unrecognized operator for disorder operator");

    psi.position(start);

    itensor::ITensor disorder = psi(start);
    disorder *= op(sites, op_type, start);
    disorder *= dag(prime_inds(psi(start), "Site", rightLinkIndex(psi, start)));

    for(int pos=start+1; pos < stop; pos++)
    {
        disorder *= psi(pos);
        disorder *= op(sites, op_type, pos);
        disorder *= dag(prime_inds(psi(pos), "Site", "Link"));
    }

    disorder *= psi(stop);
    disorder *= op(sites, op_type, stop);
    disorder *= dag(prime_inds(psi(stop), "Site", leftLinkIndex(psi, stop)));

    return disorder;
}

template<unsigned int N>
Real
clocks::compute_disorder(const clocks::Clock<N> & sites, itensor::MPS psi, const string & op_type, const Interval & interv)
{
    return elt(compute_disorder_IT(sites, psi, op_type, interv));
}

template<unsigned int N>
Complex
clocks::compute_disorderC(const clocks::Clock<N> & sites, itensor::MPS psi, const string & op_type, const Interval & interv)
{
    /* auto a = eltC(compute_disorder_IT(sites, psi, op_type, interv)); */
    /* print(a); */
    /* return a; */
    return eltC(compute_disorder_IT(sites, psi, op_type, interv));
}

#endif
