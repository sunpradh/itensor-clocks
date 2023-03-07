#ifndef __CLOCK_CLASS_H
#define __CLOCK_CLASS_H

#include <cmath>
#include "itensor/all.h"
#include "types.h"

namespace clocks {
///
/// Z_N clock models (N template parameter)
///
template<unsigned N>
class ClockSite {
private:
    it::Index s;
    complex omega = exp(complex(0, 2 * M_PI / N));

public:
    ClockSite(const it::Index & I) : s(I) {};
    ClockSite(const it::Args & args = it::Args::global());

    it::Index index() const { return s; }
    it::IndexVal state(const string & state) const;

    // Operators: "X", "Xdag", "Z", "Zdag"
    it::ITensor op(
        const string & opname,
        const it::Args & args = it::Args::global()
    ) const;
};

// Contains also a specialization for N=2, where X and Z are Hermitian

// SiteSet subclass with ClockSite type
template<unsigned int N>
using Clock = it::BasicSiteSet<ClockSite<N>>;

// Check if string represent a valid op
inline bool is_valid_op(const string & op){
    return op == "X" || op == "Xdag" || op == "Z" || op == "Zdag";
}

/************************************************************/

///
/// Default constructor for the ClockSite
///
template<unsigned N>
ClockSite<N>::ClockSite(const it::Args & args)
{
    auto ts = it::TagSet(it::format("Site,Z%d", N));
    if (args.defined("SiteNumber"))
        ts.addTags("n=" + it::str(args.getInt("SiteNumber")));
    if(args.getBool("ConserveQNs", true)) {
        // TODO gotta find a better way
        switch(N) {
            case 3:
                s = it::Index(
                        it::QN({"T", 0, 3}), 1,
                        it::QN({"T", 1, 3}), 1,
                        it::QN({"T", 2, 3}), 1,
                        it::Out, ts
                    );
                break;
            case 4:
                s = it::Index(
                        it::QN({"T", 0, 4}), 1,
                        it::QN({"T", 1, 4}), 1,
                        it::QN({"T", 2, 4}), 1,
                        it::QN({"T", 3, 4}), 1,
                        it::Out, ts
                    );
                break;
            case 5:
                s = it::Index(
                        it::QN({"T", 0, 5}), 1,
                        it::QN({"T", 1, 5}), 1,
                        it::QN({"T", 2, 5}), 1,
                        it::QN({"T", 3, 5}), 1,
                        it::QN({"T", 4, 5}), 1,
                        it::Out, ts
                    );
                break;
            case 6:
                s = it::Index(
                        it::QN({"T", 0, 6}), 1,
                        it::QN({"T", 1, 6}), 1,
                        it::QN({"T", 2, 6}), 1,
                        it::QN({"T", 3, 6}), 1,
                        it::QN({"T", 4, 6}), 1,
                        it::QN({"T", 5, 6}), 1,
                        it::Out, ts
                    );
                break;
            default:
                throw it::ITError("Case N > 6 not supported, modify clock.h if you need.");
                break;
        }
    } else {
    s = it::Index(N,ts);
}
}

/// Get the state of the clock site
template<unsigned N>
it::IndexVal ClockSite<N>::state(
    const string & state
) const {
    // state name is "n" for n=0,...,N-1
    unsigned n = std::stoi(state);
    if (n >= 0 && n < N)
        return s(n+1);
    throw it::ITError("State " + state + " not recognized");
    return it::IndexVal{};
}

template<unsigned N>
it::ITensor ClockSite<N>::op(
    const string   & opname,
    const it::Args & args
) const {
    // auto s = this->s;
    auto sP = it::prime(s);
    auto Op = it::ITensor(it::dag(s), sP);
    // we can do better
    if (opname == "Z") {
        for (auto i : it::range1(N))
        Op.set(s(mod1(i+1, N)), sP(i), 1.0);
    } else if(opname == "Zdag") {
        for (auto i : it::range1(N))
        Op.set(s(i), sP(mod1(i+1, N)), 1.0);
    } else if(opname == "X") {
        for (auto i : it::range1(N))
        Op.set(s(i), sP(i), pow(omega, i-1));
    } else if(opname == "Xdag") {
        for (auto i : it::range1(N))
        Op.set(s(i), sP(i), pow(omega, N-i+1));
    } else {
        throw it::ITError("Operator \"" + opname + "\" name not recognized");
    }
    return Op;
}

template<>
class ClockSite<2> {
    it::Index s;
public:

    ClockSite(const it::Index & I) : s(I) { }
    ClockSite(const it::Args & args = it::Args::global()) {
        auto ts = it::TagSet("Site,Z2");
        if (args.defined("SiteNumber"))
            ts.addTags("n=" + it::str(args.getInt("SiteNumber")));
        s = it::Index(2,ts);
    }

    it::Index index() const { return s; }

    it::IndexVal state(const string & state) const {
        if (state == "Up") {
            return s(1);
        } else if (state == "Dn") {
            return s(2);
        } else {
            throw it::ITError("State " + state + " not recognized");
        }
        return it::IndexVal{};
    }

    it::ITensor op(const string & opname, const it::Args & args = it::Args::global()) const {
        auto sP = it::prime(s);
        auto Op = it::ITensor(it::dag(s), sP);
        if (opname == "X" or opname == "Xdag") {
            Op.set(s(1), sP(2), +1.0);
            Op.set(s(2), sP(1), +1.0);
        } else if (opname == "Z" or opname == "Zdag") {
            Op.set(s(1), sP(1), +1.0);
            Op.set(s(2), sP(2), -1.0);
        } else {
            throw it::ITError("Operator \"" + opname + "\" name not recognized");
        }
        return Op;
    }

};

}


#endif
