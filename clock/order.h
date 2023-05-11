#ifndef __CLOCK_ORDER_H
#define __CLOCK_ORDER_H

#include <numeric>
#include <vector>

#include "itensor/all.h"
#include "clock.h"
#include "../utils/ranges.h"

template<typename T> using vector = std::vector<T>;

/************************************************************/
namespace clocks {

namespace ut = utils;

/// Computer order parameter
/// yields real result
template<unsigned N, typename... Args>
double
compute_order(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
);

/// Computer order parameter
/// yields complex result
template<unsigned N, typename... Args>
complex
compute_orderC(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
);

/// Computer order parameter
/// yields real result
template<unsigned N, typename... Args>
double
compute_bulk_order(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
);

/// Computer order parameter
/// yields complex result
template<unsigned N, typename... Args>
complex
compute_bulk_orderC(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
);

/************************************************************/

template<unsigned N, typename... Args>
double compute_order(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
) {
    auto ops = vector<string>{op_types...};
    vector<double> results;
    results.reserve(ops.size());
    for (const auto & expt : it::expect(psi, sites, ops))
        results.push_back(
            std::accumulate(expt.begin(), expt.end(), .0)
        );
    return std::accumulate(results.begin(), results.end(), .0) / double(it::length(sites));
}

template<unsigned N, typename... Args>
complex compute_orderC(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
) {
    auto ops = vector<string>{op_types...};
    vector<complex> results;
    results.reserve(ops.size());
    for (const auto & expt : it::expectC(psi, sites, ops))
        results.push_back(
            std::accumulate(expt.begin(), expt.end(), complex(.0))
        );
    return std::accumulate(results.begin(), results.end(), complex(.0)) / double(it::length(sites));
}

// Bulk order

template<unsigned N, typename... Args>
double compute_bulk_order(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
) {
    auto ops = vector<string>{op_types...};
    vector<double> results;
    results.reserve(ops.size());

    int L = it::length(sites);
    auto site_list = ut::range(L/4, (3*L/4) + 1).to_vector();

    for (const auto & expt : it::expect(psi, sites, ops, site_list))
        results.push_back(
            std::accumulate(expt.begin(), expt.end(), .0)
        );
    return std::accumulate(results.begin(), results.end(), .0) / double(site_list.size());
}

template<unsigned N, typename... Args>
complex compute_bulk_orderC(
    const Clock<N> & sites,
    it::MPS & psi,
    Args... op_types
) {
    auto ops = vector<string>{op_types...};
    vector<complex> results;
    results.reserve(ops.size());

    int L = it::length(sites);
    auto site_list = ut::range(L/4, (3*L/4) + 1).to_vector();

    for (const auto & expt : it::expectC(psi, sites, ops, site_list))
        results.push_back(
            std::accumulate(expt.begin(), expt.end(), complex(.0))
        );
    return std::accumulate(results.begin(), results.end(), complex(.0)) / double(site_list.size());
}

}
#endif
