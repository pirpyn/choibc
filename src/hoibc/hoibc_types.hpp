#ifndef _H_HOIBC_TYPES
#define _H_HOIBC_TYPES

#include <complex>
#include <array>
#include <valarray>

#if (__cplusplus >= 202002L )
#define _HOIBC_HAS_CPP20
#endif
#if (__cplusplus >= 201703L )
#define _HOIBC_HAS_CPP17
#endif

namespace hoibc
{
  using integer = long int;
  using real    = double;
  using complex = std::complex<real>;

  template<typename T> using array = std::valarray<T>;
  template<typename T> using matrix = std::array<std::array<T,2>,2>;
  template<typename T> using big_matrix = array<array<matrix<T>>>;

  static big_matrix<real> empty_bigmatrix_real;

}

#endif