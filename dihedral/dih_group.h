#ifndef __DIHEDRAL_GROUP_CLASS_H
#define __DIHEDRAL_GROUP_CLASS_H

#include <array>
#include <vector>
#include <utility>
#include "dihedral.h"

using uint = unsigned int;
using uint_pair = std::pair<uint, uint>;


template<uint N>
class DihedralGroup
{
    private:
    // store the elements of DN
    typedef Dihedral<N> element_t;
    using ElementsArray = std::array<Dihedral<N>, 2*N>;
    ElementsArray elements;

    public:
    DihedralGroup()
    {
        // generate all the elements
        for (int s=0; s < 2; s++)
            for (int r=0; r < N; r++)
                elements[r + s*N] = Dihedral<N>(r, s);
    }

    // iterators
    auto begin()  { return elements.begin(); }
    auto end()    { return elements.end(); }
    auto cbegin() const { return elements.cbegin(); }
    auto cend()   const { return elements.cend(); }
    auto size()   { return 2*N; }

};

template<uint N>
class DihedralGroupIrrep
{
    // TODO

};

#endif
