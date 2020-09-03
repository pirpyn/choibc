#ifndef _H_HOIBC_TYPES
#define _H_HOIBC_TYPES

#include <complex>
#include <array>
#include <valarray>

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

/* To remove when hitting v1.0 */
#define NOTFINISHED(where) std::cerr << where << ": not finished." << std::endl;

#endif