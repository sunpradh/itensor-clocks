#ifndef __CLOCK_ORDER_H
#define __CLOCK_ORDER_H

#include "clock.h"
#include <functional>
#include <numeric>

using std::string;
using Real = itensor::Real;
// using itensor::ITensor;
using itensor::MPS;
using itensor::MPO;
using itensor::range1;

/************************************************************/
namespace clocks {

template<unsigned N>
MPO orderMPO(
    const Clock<N> & sites,
    const string   & op_type
);

template<unsigned N>
Real compute_order(
    const Clock<N> & sites,
    const MPS      & psi,
    const string   & op_type
);

template<unsigned N>
Complex compute_orderC(
    const Clock<N> & sites,
    const MPS      & psi,
    const string   & op_type
);

}
/************************************************************/


template<unsigned N>
MPO
clocks::orderMPO(
    const clocks::Clock<N> & sites,
    const string & op_type
)
{
    int L = length(sites);
    auto order_ampo = itensor::AutoMPO(sites);
    if (!clocks::is_valid_op(op_type))
        throw std::runtime_error("Unrecognized operator type for order operator");

    if (op_type == "X" or op_type == "Xdag")
        for (auto i : range1(L))
        {
            order_ampo += 0.5, "X", i;
            order_ampo += 0.5, "Xdag", i;
        }

    if (op_type == "Z" or op_type == "Zdag")
        for (auto i : range1(L))
        {
            order_ampo += 0.5, "Z", i;
            order_ampo += 0.5, "Zdag", i;
        }

    auto order_mpo = itensor::toMPO(order_ampo);
    order_mpo /= L;
    return order_mpo;
}

template<unsigned N>
Real
clocks::compute_order(
    const clocks::Clock<N> & sites,
    const MPS & psi,
    const string & op_type
)
{
    return inner(psi, orderMPO(sites, op_type), psi);
}

template<unsigned N>
Complex
clocks::compute_orderC(
    const clocks::Clock<N> & sites,
    const MPS & psi,
    const string & op_type
)
{
    return innerC(psi, orderMPO(sites, op_type), psi);
}

#endif
