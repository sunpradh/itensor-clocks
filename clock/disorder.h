#ifndef __CLOCK_DISORDER_H
#define __CLOCK_DISORDER_H

#include "clock.h"
#include "types.h"
#include "../utils/all.h"

/************************************************************/

namespace clocks {

template<unsigned int N>
it::ITensor
compute_disorder_IT(
    const Clock<N> & sites,
    it::MPS   &    psi,
    const string   & op_type,
    const utils::Interval & interv
);


template<unsigned int N>
double
compute_disorder(
    const Clock<N> & sites,
    it::MPS   &    psi,
    const string   & op_type,
    const utils::Interval & interv
);

template<unsigned int N>
complex
compute_disorderC(
    const Clock<N> & sites,
    it::MPS   &    psi,
    const string   & op_type,
    const utils::Interval & interv
);

/************************************************************/

// Compute disorder parameter
// Return a scalar ITensor
template<unsigned int N>
it::ITensor compute_disorder_IT(
    const Clock<N> & sites,
    it::MPS & psi,
    const string & op_type,
    const utils::Interval & interv
) {
    if (!is_valid_op(op_type))
        throw std::runtime_error("Unrecognized operator for disorder operator");

    int L = length(sites);
    auto [begin, end] = interv;
    if (begin < 0 || end > L || begin >= end)
        throw std::runtime_error("Incorrect range for disorder operator");

    psi.position(begin);

    it::ITensor disorder = psi(begin);
    disorder *= op(sites, op_type, begin);
    disorder *= dag(utils::prime_inds(psi(begin), "Site", it::rightLinkIndex(psi, begin)));

    for(int pos=begin+1; pos < end; pos++) {
        disorder *= psi(pos);
        disorder *= op(sites, op_type, pos);
        disorder *= dag(utils::prime_inds(psi(pos), "Site", "Link"));
    }

    disorder *= psi(end);
    disorder *= op(sites, op_type, end);
    disorder *= dag(utils::prime_inds(psi(end), "Site", it::leftLinkIndex(psi, end)));

    return disorder;
}

template<unsigned int N>
double compute_disorder(
    const Clock<N> & sites,
    it::MPS & psi,
    const string & op_type,
    const utils::Interval & interv
) {
    return elt(compute_disorder_IT(sites, psi, op_type, interv));
}

template<unsigned int N>
complex compute_disorderC(
    const Clock<N> & sites,
    it::MPS & psi,
    const string & op_type,
    const utils::Interval & interv
) {
    return eltC(compute_disorder_IT(sites, psi, op_type, interv));
}

}
#endif
