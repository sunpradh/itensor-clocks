#ifndef __CLOCK_UTILS_TIME
#define __CLOCK_UTILS_TIME

#include <chrono>
#include <string>
#include <type_traits>

namespace utils {

namespace chrono = std::chrono;

namespace time {
    using ns = chrono::nanoseconds;
    using us = chrono::microseconds;
    using ms = chrono::milliseconds;
    using s = chrono::seconds;
    using m = chrono::minutes;
    using h = chrono::hours;
}

// returns the unit of measure string for Time_t
template<typename Time_t>
constexpr auto units() {
    if constexpr (std::is_same<Time_t, time::ns>::value) return "ns";
    if constexpr (std::is_same<Time_t, time::us>::value) return "μs";
    if constexpr (std::is_same<Time_t, time::ms>::value) return "ms";
    if constexpr (std::is_same<Time_t, time::s>::value)  return "s";
    if constexpr (std::is_same<Time_t, time::m>::value)  return "m";
    if constexpr (std::is_same<Time_t, time::h>::value)  return "h";
    return "";
}

}

#endif
