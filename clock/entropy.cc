#include "entropy.h"

using Real = itensor::Real;
Real entropy_vN(itensor::MPS & psi, int pos, const itensor::Args & args)
{
    psi.position(pos);
    auto link_ind = leftLinkIndex(psi, pos);
    auto site_ind = siteIndex(psi, pos);

    auto [U,S,V] = svd(psi(pos), {link_ind, site_ind});
    auto common_ind = commonIndex(U,S);

    auto cutoff = args.getReal("Cutoff", 1E-12);
    Real SvN = 0;
    for(auto n : itensor::range1(dim(common_ind)))
    {
        auto p = sqrt(elt(S, n, n));
        if (p > cutoff) SvN += -p*log(p);
    }

    return SvN;
}
