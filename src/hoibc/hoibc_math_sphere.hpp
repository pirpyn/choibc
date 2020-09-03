#ifndef _H_HOIBC_MATH_SPHERE
#define _H_HOIBC_MATH_SPHERE

#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc {
  namespace sphere {
    big_matrix<complex> impedance_infinite(const array<real>& n, const real& k0, const material_t& material, const real& inner_radius);

    big_matrix<complex> reflexion_infinite(const array<real>& n, const real& k0, const material_t& material, const real& inner_radius);
  }
}

#endif