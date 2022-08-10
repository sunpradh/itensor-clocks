#ifndef __CLOCK_UTILS_H
#define __CLOCK_UTILS_H

#include <vector>
#include <cmath>
#include <string>

#include "itensor/mps/mps.h"
#include "itensor/indexset.h"

/************************************************************/

using string = std::string;
using Interval = std::pair<int, int>;
using size_t = std::size_t;

constexpr int mod1(int i, int n) { return ((i - 1) % n) + 1; }

// Shorthand for ranges and linspaces like in numpy
#include "ranges.h"

// dumping contents
#include "dump.h"

// Benchmarks
#include "benchmark.h"

// Shorthand for priming multiple indices at once
template<typename T, typename... Args>
itensor::ITensor prime_inds(itensor::ITensor psi, T ind, Args... args);

template<typename T>
itensor::ITensor prime_inds(itensor::ITensor psi, T ind);

// progress bar
class ProgressBar;

// Converto unordered_map to vector
template<typename T1, typename T2>
std::vector<T2>
umap_to_vector(std::unordered_map<T1, T2> in, std::vector<T1> keys)
{
    std::vector<T2> out;
    for (auto k : keys)
    {
        out.push_back(in[k]);
    }
    return out;
}

/************************************************************/


template<typename T>
itensor::ITensor prime_inds(itensor::ITensor psi, T ind)
{
    return itensor::prime(psi, ind);
}

template<typename T, typename... Args>
itensor::ITensor prime_inds(itensor::ITensor psi, T ind, Args... args)
{
    return prime_inds(itensor::prime(psi, ind), args...);
}

class ProgressBar {
    private:
    unsigned int counter = 0;
    unsigned int total;

    public:
    ProgressBar() : total(0) {};
    ProgressBar(unsigned int total_) : total(total_) {};
    ProgressBar(unsigned int start, unsigned int total_) : counter(start), total(total_) {};

    void increment() {
        if (counter < total)
            counter++;
    }

    void print() {
        printf("\r%d/%d", counter, total);
    }

};


#endif
