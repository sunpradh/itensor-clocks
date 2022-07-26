#ifndef __CLOCK_ORDER_H
#define __CLOCK_ORDER_H

#include "clock.h"

using Real = itensor::Real;
/************************************************************/
namespace clocks {

template<unsigned int N>
itensor::MPO orderMPO(
    const clocks::Clock<N> & sites,
    const std::string      & type
);

template<unsigned int N>
Real compute_order(
    const clocks::Clock<N> & sites,
    const itensor::MPS     & psi,
    const std::string      & type
);

template<unsigned int N>
Complex compute_orderC(
    const clocks::Clock<N> & sites,
    const itensor::MPS     & psi,
    const std::string      & type
);

}
/************************************************************/


template<unsigned int N>
itensor::MPO
clocks::orderMPO(const clocks::Clock<N> & sites, const std::string & type)
{
    int L = length(sites);
    auto order_ampo = itensor::AutoMPO(sites);
    if (type == "X")
        for (int i=1; i<=L; i++)
        {
            order_ampo += 0.5, "X", i;
            order_ampo += 0.5, "Xdag", i;
        }
    else
        if (type == "Z")
            for (int i=1; i<=L; i++)
            {
                order_ampo += 0.5, "Z", i;
                order_ampo += 0.5, "Zdag", i;
            }
        else
            throw std::runtime_error("Unrecognized type for OrderOperator");

    auto order_mpo = toMPO(order_ampo);
    order_mpo /= L;
    return order_mpo;
}

template<unsigned int N>
Real
clocks::compute_order(const clocks::Clock<N> & sites, const itensor::MPS & psi, const std::string & type)
{
    return inner(psi, orderMPO(sites, type), psi);
}

template<unsigned int N>
Complex
clocks::compute_orderC(const clocks::Clock<N> & sites, const itensor::MPS & psi, const std::string & type)
{
    return innerC(psi, orderMPO(sites, type), psi);
}

#endif
