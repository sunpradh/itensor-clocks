#ifndef __DIHEDRAL_CLASS_H
#define __DIHEDRAL_CLASS_H

#include <array>
#include <utility>
#include <string>


using uint = unsigned int;
using uint_pair = std::pair<uint, uint>;

template<uint N>
class Dihedral
{
    private:
    /*
    * Any element of the Dihedral group can be written as
    *       r**k s**j,      with k=0,...,N-1 and j=0,1
    * Obviously:
    *       r**N = e (unit element)
    *       s**2 = e
    *       s r s = r**(-1)    or    r s = s r**(-1)
    */
    uint r;
    uint s;

    public:
    // default constructor: unit element
    Dihedral() : r(0), s(0) {};
    // basic constructor
    Dihedral(uint r_, uint s_) : r(r_ % N), s(s_ % 2) {};
    // pair constructor
    Dihedral(uint_pair p) : r(p.first % N), s(p.second % 2) {};
    // copy construct
    Dihedral(const Dihedral & h) : r(h.r % N), s(h.s % 2) {};
    // assignment
    Dihedral operator = (const Dihedral & h) {
        this->r = h.r;
        this->s = h.s;
        return *this;
    }

    // check if is the unit element
    bool is_unit() const
    {
        return (r == 0 && s == 0);
    }

    // group inversion
    Dihedral inv() const {
        return s==0 ? Dihedral(N-r, 0) : (*this);
    }

    // Group multiplication
    Dihedral operator * (const Dihedral & h) const
    {
        // multiplication by unit is trivial
        /* if (h.isunit())     return (*this); */
        /* if (this->isunit()) return h; */

        if (s==1)
            return Dihedral(r + (N - h.r), (s + h.s) % 2);
        else
            return Dihedral(r + h.r, h.s);

    }

    // Comparison operators
    bool operator == (const Dihedral & h) const
    {
        return (r == h.r && s == h.s);
    }
    bool operator != (const Dihedral & h) const
    {
        return (r != h.r || s != h.s);
    }

    // convert to string
    std::string str() const {
        if (this->is_unit())
            return "e";
        else if  (r == 0 && s == 1)
            return "s";
        else if (s == 0 && r != 0)
            return r == 1 ? "r" : "r" + std::to_string(r);
        else
            return (r == 1 ? "r" : "r" + std::to_string(r)) + "s";
    }

};

template<uint N>
bool operator == (const std::string & str, const Dihedral<N> g)
{
    return str == g.str();
}

template<uint N>
bool operator == (const Dihedral<N> g, const std::string & str)
{
    return str == g.str();
}


#endif
