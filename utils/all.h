#ifndef __CLOCK_UTILS_H
#define __CLOCK_UTILS_H

#include <string>
#include <vector>
#include <unordered_map>

#include "itensor/all.h"


/************************************************************/
// List of useful functions, that do not directly depends
// on the Clock<N> class
/************************************************************/

namespace utils {

namespace it = itensor;
using string = std::string;
using Interval = std::pair<int, int>;

// modulus counting from 1
//      mod1(N, N) = N,
//      mod1(N + 1, N) = 1,
// different from the usual mod operation
constexpr int mod1(int i, int n) {
    return ((i - 1) % n) + 1;
}

// Primes one index
template<typename T>
it::ITensor prime_inds(it::ITensor psi, T ind) {
    return it::prime(psi, ind);
}

// Shorthand for priming multiple indices at once
template<typename T, typename... Args>
it::ITensor
prime_inds(it::ITensor psi, T ind, Args... args) {
    return prime_inds(it::prime(psi, ind), args...);
}

// Convert unordered_map to vector
template<typename T1, typename T2>
std::vector<T2>
umap_to_vector(std::unordered_map<T1, T2> in, std::vector<T1> keys) {
    std::vector<T2> out;
    for (auto k : keys)
        out.push_back(in[k]);
    return out;
}

}

// Shorthand for ranges and linspaces like in numpy
#include "ranges.h"

// Table class
#include "table.h"

// Benchmarks
#include "benchmark.h"

// Timers
#include "timer.h"

/************************************************************/


#endif
