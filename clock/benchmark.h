#include <iostream>
#include <chrono>
#include <functional>
#include <numeric>
#include <utility>
#include <type_traits>

using namespace std::chrono;
using func_t = std::function<void()>;


template<typename Time_t>
constexpr std::string
units_suffix() {
    if (std::is_same<Time_t, nanoseconds>::value) return "ns";
    if (std::is_same<Time_t, microseconds>::value) return "us";
    if (std::is_same<Time_t, milliseconds>::value) return "ms";
    if (std::is_same<Time_t, seconds>::value) return "s";
}


template<unsigned N_iter = 1>
class benchmark {
    private:
    std::array<time_point<high_resolution_clock>, N_iter+1> time_points;
    func_t func;

    public:
    benchmark(const func_t & func_) {
        func = func_;
        time_points.front() = high_resolution_clock::now();
        for (unsigned i=0; i<N_iter; i++) {
            func();
            time_points.at(i+1) = high_resolution_clock::now();
        }
    };

    // total duration
    template<typename Time_t, typename Rep = long>
    auto total_duration() {
        return duration_cast<Time_t, Rep>(time_points.back() - time_points.front()).count();
    };

    // single lap duration
    template<typename Time_t, typename Rep = long>
    auto duration(unsigned i) {
        return duration_cast<Time_t, Rep>(time_points.at(i+1) - time_points.at(i)).count();
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
        auto unit = units_suffix<Time_t>();
        std::cout
            << "Number of iterations: " << N_iter << "\n"
            << "Total duration: " << total << " " << unit << "\n"
            << "Average duration: (" << avg << " ± " << dev << ") " << unit << "\n";
    }
};
