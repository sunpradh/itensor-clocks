#ifndef __CLOCK_CLASS_H
#define __CLOCK_CLASS_H

#include <cmath>
#include "itensor/all.h"
#include "utils.h"

using std::string;
using itensor::Args;
using itensor::Complex;
using itensor::ITensor;
using itensor::Index;
using itensor::IndexVal;
using itensor::range1;
using itensor::QN;

/************************************************************/

// Z_N clock models (N template parameter)
namespace clocks {
// Operators: "X", "Xdag", "Z", "Zdag"
template<unsigned N>
class ClockSite {
private:
    Index s;
    Complex omega = exp(Complex(0, 2 * M_PI / N));

public:
    ClockSite(Index const & I) : s(I) {};
    ClockSite(Args const & args = Args::global());

    Index index() const { return s; }
    IndexVal state(string const & state);

    ITensor op(
        string const & opname,
        Args const & args = Args::global()
    ) const;
    static bool is_valid_op(string const & op);
};

// Specialization for N=2, where X and Z are Hermitian
// template<> class ClockSite<2>;

template<unsigned int N>
using Clock = itensor::BasicSiteSet<ClockSite<N>>;
}
/************************************************************/

template<unsigned N>
clocks::ClockSite<N>::ClockSite(Args const & args)
{
    auto ts = itensor::TagSet(itensor::format("Site,Z%d", N));
    if (args.defined("SiteNumber"))
        ts.addTags("n=" + itensor::str(args.getInt("SiteNumber")));
    if(args.getBool("ConserveQNs", true)) {
        // gotta find a better way
        switch(N) {
            case 3:
                s = Index(
                        QN({"T", 0, 3}), 1,
                        QN({"T", 1, 3}), 1,
                        QN({"T", 2, 3}), 1,
                        itensor::Out, ts
                    );
                break;
            case 4:
                s = Index(
                        QN({"T", 0, 4}), 1,
                        QN({"T", 1, 4}), 1,
                        QN({"T", 2, 4}), 1,
                        QN({"T", 3, 4}), 1,
                        itensor::Out, ts
                    );
                break;
            case 5:
                s = Index(
                        QN({"T", 0, 5}), 1,
                        QN({"T", 1, 5}), 1,
                        QN({"T", 2, 5}), 1,
                        QN({"T", 3, 5}), 1,
                        QN({"T", 4, 5}), 1,
                        itensor::Out, ts
                    );
                break;
            case 6:
                s = Index(
                        QN({"T", 0, 6}), 1,
                        QN({"T", 1, 6}), 1,
                        QN({"T", 2, 6}), 1,
                        QN({"T", 3, 6}), 1,
                        QN({"T", 4, 6}), 1,
                        QN({"T", 5, 6}), 1,
                        itensor::Out, ts
                    );
                break;
            default:
                throw itensor::ITError("Case N > 6 not supported, modify clock.h if you need.");
                break;
            }
    } else {
        s = itensor::Index(N,ts);
    }
}

template<unsigned N>
itensor::IndexVal
clocks::ClockSite<N>::state(
    string const & state
)
{
    // state name is "n" for n=0,...,N-1
    unsigned n = std::stoi(state);
    if (n >= 0 && n < N)
        return s(n+1);
    throw itensor::ITError("State " + state + " not recognized");
    return IndexVal{};
}

template<unsigned N>
itensor::ITensor
clocks::ClockSite<N>::op(
    string const & opname,
    Args   const & args
) const
{
    // auto s = this->s;
    auto sP = itensor::prime(s);
    auto Op = ITensor(itensor::dag(s), sP);
    // we can do better
    if (opname == "Z") {
        for (auto i : range1(N))
            Op.set(s(mod1(i+1, N)), sP(i), 1.0);
    } else if(opname == "Zdag") {
        for (auto i : range1(N))
            Op.set(s(i), sP(mod1(i+1, N)), 1.0);
    } else if(opname == "X") {
        for (auto i : range1(N))
            Op.set(s(i), sP(i), pow(omega, i-1));
    } else if(opname == "Xdag") {
        for (auto i : range1(N))
            Op.set(s(i), sP(i), pow(omega, -i+1));
    } else {
        throw itensor::ITError("Operator \"" + opname + "\" name not recognized");
    }
    return Op;
}

template<unsigned N>
bool
clocks::ClockSite<N>::is_valid_op(
    string const & op
)
{
    return op == "X" || op == "Xdag" || op == "Z" || op == "Zdag";
}


template<>
class clocks::ClockSite<2>
{
    Index s;
    public:

    ClockSite(Index const& I) : s(I) { }
    ClockSite(Args const& args = itensor::Args::global())
    {
        auto ts = itensor::TagSet("Site,Z2");
        if (args.defined("SiteNumber"))
            ts.addTags("n=" + itensor::str(args.getInt("SiteNumber")));
        s = Index(2,ts);
    }

    Index index() const { return s; }

    IndexVal state(string const& state) {
        if (state == "Up") {
            return s(1);
        } else if (state == "Dn") {
            return s(2);
        } else {
            throw itensor::ITError("State " + state + " not recognized");
        }
        return IndexVal{};
    }

    ITensor op(string const& opname, itensor::Args const& args = itensor::Args::global()) const
    {
        auto sP = itensor::prime(s);
        auto Op = ITensor(itensor::dag(s), sP);
        if (opname == "X" or opname == "Xdag") {
            Op.set(s(1), sP(2), +1.0);
            Op.set(s(2), sP(1), +1.0);
        } else if (opname == "Z" or opname == "Zdag") {
            Op.set(s(1), sP(1), +1.0);
            Op.set(s(2), sP(2), -1.0);
        } else {
            throw itensor::ITError("Operator \"" + opname + "\" name not recognized");
        }
        return Op;
    }

    bool is_valid_op(string const & op)
    {
        return op == "X" || op == "Xdag" || op == "Z" || op == "Zdag";
    }


};


#endif
