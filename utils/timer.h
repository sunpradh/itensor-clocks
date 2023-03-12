#ifndef __CLOCK_UTILS_TIMER
#define __CLOCK_UTILS_TIMER

#include <chrono>
#include <iostream>
#include <string>

#include "time.h"

namespace utils {

class Timer {
    using clock = chrono::steady_clock;
    using time_point = chrono::time_point<clock>;
    time_point start_{}, stop_{};

public:
    Timer reset() {
        start_ = time_point{};
        stop_ = time_point{};
        return *this;
    }

    Timer start() {
        start_ = clock::now();
        return *this;
    }

    Timer stop() {
        stop_ = clock::now();
        return *this;
    }

    template<typename Time_t, typename Rep = long>
    Rep duration() const {
        if (start_ == time_point{} || start_ > stop_)
            throw "Timer not started";
        if (stop_ == time_point{})
            throw "Timer not stopped";
        return chrono::duration_cast<Time_t, Rep>(stop_ - start_).count();
    }

    template<typename Time_t, typename Rep = long>
    std::string str() {
        auto dur = this->duration<Time_t, Rep>();
        return std::to_string(dur)
                + std::string(" ")
                + std::string(units<Time_t>());
    }

};


std::ostream &
operator<<(std::ostream & output, const Timer & timer) {
    auto time_ms = timer.duration<time::s>();
    if (time_ms == 0) {
        output << "<1ms";
        return output;
    }

    auto time_s = timer.duration<time::s>();
    if (time_s == 0) {
        output << time_ms << units<time::ms>();
        return output;
    }

    auto time_h = timer.duration<time::h>();
    auto time_m = timer.duration<time::m>();
    if (time_h > 0)
        output << time_h << units<time::h>();
    if (time_m > 0)
        output << time_m << units<time::m>();
    output << time_s << units<time::s>();
    return output;

}

}


#endif
