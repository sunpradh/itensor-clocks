#ifndef __CLOCK_UTILS_RANGES_H
#define __CLOCK_UTILS_RANGES_H

#include <cmath>
#include <vector>
#include <array>
#include <stdexcept>


/************************************************************/
namespace utils {

template<typename T> class range;

template<typename T> class linspace;

/************************************************************/

template<typename T>
class range_iter {
protected:
    T current, step;

public:
    range_iter(T current_, T step_ = T(1)) : current(current_), step(step_) {};
    range_iter() : current(T(0)), step(T(1)) {};

    T operator*() const { return current; }
    T const* operator->() const { return current; }

    range_iter& operator++() { current += step; return *this; }

    bool operator==(const range_iter & other) {
        // consider also cases where the step is non-integer or negative (i think)
        if (step > 0)
            return current >= other.current;
        else
            return current < other.current;
    }
};


template<typename T>
class range {
public:
    using iterator = range_iter<T>;

    range(T begin, T end, T step = T(1)) : begin_(begin, step), end_(end, step) {}
    range(T end) : range(0, end, 1) {}

    iterator begin() { return begin_; }
    iterator end() { return end_; }

private:
    iterator begin_, end_;
};


template<typename T>
class linspace {
public:
    class iterator : public range_iter<T> {
    public:
        iterator(T current_, T step_ = T(1)) : range_iter<T>(current_, step_) {}
        iterator() : range_iter<T>(T(0), T(1)) {}

        bool operator==(const iterator & other) {
            if (this->step > 0)
                return this->current > other.current;
            else
                return this->current < other.current;
        }
    };

    linspace(T begin, T end, std::size_t npoints) : npoints_(npoints) {
        auto step = (end - begin) / T(npoints-1);
        begin_ = iterator(begin, T(step));
        end_ = iterator(end, T(step));
    }

    iterator begin() { return begin_; }
    iterator end() { return end_; }
    std::size_t size() { return npoints_; }

    std::vector<T> to_vector() {
        std::vector<T> vec;
        for (auto x : *this) vec.push_back(x);
        return vec;
    }

    template<size_t N>
    std::array<T, N> to_array() {
        if (N != size())
            throw std::invalid_argument("Array of the wrong size");
        std::array<T, N> arr;
        size_t n = 0;
        for (auto x : *this) { arr[n] = x; n++; }
        return arr;
    }


private:
    iterator begin_, end_;
    std::size_t npoints_;

};

}

#endif
