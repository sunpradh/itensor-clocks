/*
*
* DOES NOT WORK
*
*
*/
#ifndef __CLOCK_UTILS_DUMP_H
#define __CLOCK_UTILS_DUMP_H

#include <unordered_map>
#include <vector>
#include <string>
#include <ostream>
#include <iomanip>
#include <utility>

using uint    = unsigned int;
using string  = std::string;
using ostream = std::ostream;

template<typename T> using vector = std::vector<T>;
template<typename T1, typename T2> using umap = std::unordered_map<T1, T2>;
template<typename T1, typename T2> using pair = std::pair<T1, T2>;

namespace utils {

namespace dump_settings {
    // these default settings are defined in dump.cc
    extern string separator;
    extern unsigned int precision;
    extern unsigned int width;
}

template <typename Tx, typename Ty>
void
dump_single_row(
    ostream & os,
    const pair<uint, Tx> & x,
    const umap<Tx, Ty> & y
)
{
    os << std::setw(dump_settings::width) << y.at(x.second);
}

template <typename Tx, typename Ty>
void
dump_single_row(
    ostream & os,
    const pair<uint, Tx> & x,
    const vector<Ty> & y
)
{
    os << std::setw(dump_settings::width) << y[x.first];
}

template <typename Tx, typename Generic, typename ... Args>
void
dump_single_row(
    ostream & os,
    const pair<uint, Tx> & x,
    const Generic & y,
    Args ... args
)
{
    dump_single_row(os, x, y);
    if constexpr (sizeof...(args) > 1)
    {
        os << dump_settings::separator;
        dump_single_row(os, x, args...);
    }
}

template <typename Tx, typename ... Args>
void
dump(ostream & os, const vector<Tx> & x, const Args & ... args)
{
    os << std::setprecision(dump_settings::precision);
    for (uint i = 0; i < x.size(); i++)
    {
        os << std::setw(dump_settings::width) << x[i] << dump_settings::separator;
        dump_single_row(os, std::make_pair(i, x[i]), args ...);
        os << '\n';
    }
}

}

#endif
