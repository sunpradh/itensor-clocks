#ifndef __MATRIX_CLASS_H
#define __MATRIX_CLASS_H

#include <array>

template<unsigned n, typename T = double>
using Matrix = std::array<std::array<T, n>, n>;

template<unsigned n, typename T = double>
T trace(Matrix<n, T> M)
{
    auto sum = T(0);
    for (int k=0; k < n; k++) sum += M[k][k];
    return sum;
}

#endif
