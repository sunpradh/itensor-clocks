#ifndef __CLOCK_CLASS_H
#define __CLOCK_CLASS_H

#include "itensor/all.h" // to bothersome to include all the specific headers needed
#include <cmath>

using Complex = itensor::Complex;
using Real = itensor::Real;
/************************************************************/

// Z_N clock models (N template parameter)
namespace clocks {
// Operators: "X", "Xdag", "Z", "Zdag"
template<unsigned int N>
class ClockSite;

// Specialization for N=2, where X and Z are Hermitian
template<> class ClockSite<2>;

template<unsigned int N>
using Clock = itensor::BasicSiteSet<ClockSite<N>>;
}
/************************************************************/

template<unsigned int N>
class clocks::ClockSite
{
    itensor::Index s;
    Complex omega = exp(Complex(0, 2 * M_PI / N));

    public:

    ClockSite(itensor::Index const& I) : s(I) {};

    ClockSite(itensor::Args const& args = itensor::Args::global())
    {
        auto ts = itensor::TagSet(itensor::format("Site,Z%d", N));
        if (args.defined("SiteNumber"))
            ts.addTags("n=" + itensor::str(args.getInt("SiteNumber")));
        s = itensor::Index(N,ts);
    }

    itensor::Index index() const { return s; }

    itensor::IndexVal state(std::string const& state)
    {
        // state name is "n" for n=0,...,N-1
        for (auto n : itensor::range(N))
            if (state == itensor::str(n))
                return s(n+1);
        throw itensor::ITError("State " + state + " not recognized");
        return itensor::IndexVal{};
    }

    itensor::ITensor op(std::string const& opname, itensor::Args const& args = itensor::Args::global()) const
    {
        auto sP = prime(s);
        auto Op = itensor::ITensor(dag(s), sP);
        if (opname == "X")
        {
            Op.set(s(1), sP(N), 1.0);
            for (unsigned int i=1; i < N; ++i)
                Op.set(s(i+1), sP(i), 1.0);
        }
        else if(opname == "Xdag")
        {
            Op.set(s(N), sP(1), 1.0);
            for (unsigned int i=1; i < N; ++i)
                Op.set(s(i), sP(i+1), 1.0);
        }
        else if(opname == "Z")
        {
            Op.set(s(1), sP(1), 1.0);
            for (unsigned int i=2; i<=N; ++i)
                Op.set(s(i), sP(i), pow(omega, i-1));
        }
        else if(opname == "Zdag")
        {
            Op.set(s(1), sP(1), 1.0);
            for (unsigned int i=2; i<=N; ++i)
                Op.set(s(i), sP(i), pow(omega, N-i+1));
        }
        else
            throw itensor::ITError("Operator \"" + opname + "\" name not recognized");

        return Op;
    }
};

template<>
class clocks::ClockSite<2>
{
    itensor::Index s;
    public:

    ClockSite(itensor::Index const& I) : s(I) { }

    ClockSite(itensor::Args const& args = itensor::Args::global())
    {
        auto ts = itensor::TagSet("Site,Z2");
        if (args.defined("SiteNumber"))
            ts.addTags("n=" + itensor::str(args.getInt("SiteNumber")));
        s = itensor::Index(2,ts);
    }

    itensor::Index index() const { return s; }

    itensor::IndexVal state(std::string const& state)
    {
        if (state == "Up")
            return s(1);
        else
            if (state == "Dn")
                return s(2);
            else
                throw itensor::ITError("State " + state + " not recognized");
        return itensor::IndexVal{};
    }

    itensor::ITensor op(std::string const& opname, itensor::Args const& args = itensor::Args::global()) const
    {
        auto sP = prime(s);
        auto Op = itensor::ITensor(dag(s), sP);
        if (opname == "X" or opname == "Xdag")
        {
            Op.set(s(1), sP(2), +1.0);
            Op.set(s(2), sP(1), +1.0);
        }
        else
            if (opname == "Z" or opname == "Zdag")
            {
                Op.set(s(1), sP(1), +1.0);
                Op.set(s(2), sP(2), -1.0);
            }
            else
                throw itensor::ITError("Operator \"" + opname + "\" name not recognized");

        return Op;
    }

};


#endif
