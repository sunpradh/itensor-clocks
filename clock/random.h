#ifndef __CLOCK_RANDOM_H
#define __CLOCK_RANDOM_H
#include <vector>
#include <random>
#include <numeric>

#include "itensor/all.h"
#include "clock.h"

namespace clocks {

template<int Mod>
std::vector<int> random_ints_modulo(unsigned length, int sum)
{
    // Taken directly from cppreference.com
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> distrib(0, Mod-1);

    std::vector<int> vec(length);
    for (unsigned i=0; i < length-1; i++) {
        vec[i] = distrib(gen);
    }
    int partial = std::accumulate(vec.begin(), vec.end()-1, 0) % Mod;
    // Avoid negative numbers with the modulo operation
    vec.back() = (((sum - partial) % Mod) + Mod) % Mod;
    return vec;
}

template<unsigned N>
itensor::MPS
randomMPS_QN(const Clock<N> & sites, int qn)
{
    auto L = itensor::length(sites);
    auto random_ints = random_ints_modulo<N>(L, qn);
    auto state = itensor::InitState(sites);
    for (auto i : itensor::range1(L)) {
        state.set(i, std::to_string(random_ints[i-1]));
    }
    return itensor::MPS(state);
}

}


#endif
