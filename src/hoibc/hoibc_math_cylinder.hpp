#ifndef _H_HOIBC_MATH_CYLINDER
#define _H_HOIBC_MATH_CYLINDER

#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc {
  namespace cylinder {

    big_matrix<complex> impedance_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& radius);

    big_matrix<complex> reflexion_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& radius);
  }
}
#endif