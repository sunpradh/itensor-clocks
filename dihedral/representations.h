#ifndef __REPRESENTATIONS_H
#define __REPRESENTATIONS_H

#include <functional>
#include <numeric>
#include <string>
#include <unordered_map>
#include "matrix.h"

template<typename Group, typename T>
class ReprBase {
    protected:
    using element_t = typename Group::element_t;
    using function_t = std::function<T(element_t)>;
    function_t func;

    public:
    ReprBase(const function_t & func_) : func(func_) {};
    T operator () (element_t g)
    {
        return func(g);
    }

    void character(element_t);
};

template <unsigned dim, typename Group, typename T = double>
class Repr : ReprBase<Group, Matrix<dim, T> >
{
    using element_t = typename Group::element_t;

    public:
    T character(element_t g)
    {
        return trace(this->func(g));
    }
};

template <typename Group, typename T>
class Repr<1,Group,T> : ReprBase<Group, T>
{
    public:
    T character(typename Group::element_t g)
    {
        return this->func(g);
    }
};


template <unsigned dim, typename Group, typename T = double>
class ReprList
{
    using element_t = typename Group::element_t;

    std::unordered_map<std::string, Repr<dim, Group, T>> repr_list;

    public:
    // TODO


};

#endif
