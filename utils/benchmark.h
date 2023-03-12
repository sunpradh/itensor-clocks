#ifndef __CLOCK_UTILS_BENCHMARK
#define __CLOCK_UTILS_BENCHMARK

#include <iostream>
#include <chrono>
#include <functional>
#include <numeric>
#include <utility>

#include "time.h"

/************************************************************/
namespace utils {

using func_t = std::function<void()>;

// benchmark class
// it is initialized with a void function func_t, the object to benchmark
template<unsigned N_iter = 1> class Benchmark;



template<unsigned N_iter>
class Benchmark {
private:
    using clock = chrono::steady_clock;
    using time_point = chrono::time_point<clock>;
    std::array<time_point, N_iter+1> time_points;
    func_t func;

public:
    Benchmark(const func_t & func_) {
        func = func_;
        time_points.front() = clock::now();
        for (unsigned i=0; i<N_iter; i++) {
            func();
            time_points.at(i+1) = clock::now();
        }
    };

    // total duration
    template<typename Time_t, typename Rep = long>
    auto total_duration() {
        return chrono::duration_cast<Time_t, Rep>(time_points.back() - time_points.front()).count();
    };

    // single lap duration
    template<typename Time_t, typename Rep = long>
    auto duration(unsigned i) {
        return chrono::duration_cast<Time_t, Rep>(time_points.at(i+1) - time_points.at(i)).count();
    };

    // average duration
    // std deviation of durations
    template<typename Time_t, typename Rep = long>
    auto average() {
        std::array<Rep, N_iter> lap_times;
        for (unsigned i=0; i < N_iter; i++) {
            lap_times.at(i) = duration<Time_t, Rep>(i);
        }

        // Compute average
        Rep avg = std::accumulate(
                lap_times.begin(), lap_times.end(), Rep(0)
            ) / lap_times.size();

        // Compute std deviation
        Rep dev = std::accumulate(
                lap_times.begin(), lap_times.end(), Rep(0),
                [&avg](Rep sum, Rep elem){ return sum + (avg - elem)*(avg - elem); }
            );
        dev = sqrt(dev);
        return std::make_pair(avg, dev);
    };

    template<typename Time_t, typename Rep = long>
    void print_statistics() {
        auto total = total_duration<Time_t, Rep>();
        auto [avg, dev] = average<Time_t, Rep>();
        auto unit = units<Time_t>();
        std::cout
            << "\tNumber of iterations: " << N_iter << "\n"
            << "\tTotal duration: " << total << " " << unit << "\n"
            << "\tAverage duration: (" << avg << " ± " << dev << ") " << unit << "\n";
    }
};

}
#endif
