#ifndef _H_HOIBC_TYPES
#define _H_HOIBC_TYPES

#include <complex>
#include <vector>
#include <array>

namespace hoibc
{
  using integer = long int;
  using real = double;
  using complex = std::complex<real>;
  template<typename T> using matrix = std::array<std::array<T,2>,2>;
  template<typename T> using big_matrix = std::vector<std::vector<matrix<T>>>;
}
#endif