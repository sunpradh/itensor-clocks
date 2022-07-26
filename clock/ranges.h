#ifndef __CLOCK_UTILS_RANGES_H
#define __CLOCK_UTILS_RANGES_H

#include <vector>
#include <cmath>

/************************************************************/
// Shorthand for ranges and linspaces like in numpy
template<typename T>
std::vector<T> range(T end);

template<typename T>
std::vector<T> range(T begin, T end);

template<typename T>
std::vector<T> range(T begin, T end, T step);

template<typename T>
std::vector<T> linspace(T begin, T end, size_t npoints = 100);
/************************************************************/

template<typename T>
std::vector<T> range(T begin, T end, T step)
{
    std::vector<T> vec;
    vec.reserve(size_t(std::abs(end-begin)/step));
    for (T elem = begin; elem < end; elem += step)
        vec.push_back(elem);
    return vec;
}

template<typename T>
std::vector<T> range(T begin, T end)
{
    return range(begin, end, T(1));
}

template<typename T>
std::vector<T> range(T end)
{
    return range(T(0), end, T(1));
}

template<typename T>
std::vector<T> linspace(T begin, T end, size_t npoints)
{
    std::vector<T> vec;
    vec.reserve(npoints);
    T step = T(std::abs(end-begin)) / T(npoints-1);
    for (T elem = begin; elem <= end; elem += step)
        vec.push_back(elem);
    return vec;
}

#endif
