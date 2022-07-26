#ifndef __ENTROPY_H
#define __ENTROPY_H

#include "itensor/all.h"

itensor::Real entropy_vN(itensor::MPS & psi, int pos, const itensor::Args & args = itensor::Args::global());

#endif
