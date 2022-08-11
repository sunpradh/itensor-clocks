#ifndef __CLOCK_ORDER_H
#define __CLOCK_ORDER_H

#include <numeric>
#include <vector>
#include <type_traits>
#include <utility>

#include "clock.h"
#include <functional>

template<typename T> using vector = std::vector<T>;
using std::string;
using Real = itensor::Real;
using itensor::expect;
using itensor::MPS;

/************************************************************/
namespace clocks {

template<unsigned N, typename... Args>
Real
compute_order(
    const Clock<N> & sites,
    MPS & psi,
    Args... op_types
);

template<unsigned N, typename... Args>
Complex
compute_orderC(
    const Clock<N> & sites,
    MPS & psi,
    Args... op_types
);

}
/************************************************************/

template<unsigned N, typename... Args>
Real
clocks::compute_order(
    const clocks::Clock<N> & sites,
    MPS & psi,
    Args... op_types
)
{
    auto ops = vector<string>{op_types...};
    auto expts = expect(psi, sites, ops);
    vector<Real> results;
    for (const auto & expt : expts)
        results.push_back(std::accumulate(expt.begin(), expt.end(), Real(0)));
    return std::accumulate(results.begin(), results.end(), Real(0)) / double(itensor::length(sites));
}

template<unsigned N, typename... Args>
Complex
clocks::compute_orderC(
    const clocks::Clock<N> & sites,
    MPS & psi,
    Args... op_types
)
{
    auto ops = vector<string>{op_types...};
    auto expts = expectC(psi, sites, ops);
    vector<Complex> results;
    for (const auto & expt : expts)
        results.push_back(std::accumulate(expt.begin(), expt.end(), Complex(0)));
    return std::accumulate(results.begin(), results.end(), Complex(0)) / double(itensor::length(sites));
}

#endif
