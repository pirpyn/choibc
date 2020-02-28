#ifndef _H_HOIBC_TYPES
#define _H_HOIBC_TYPES

#include <complex>
#include <vector>

namespace hoibc
{
  using integer = long int;
  using real = double;
  using complex = std::complex<real>;
  template<typename T> using big_matrices = std::vector<std::vector<T[2][2]>>;
}
#endif