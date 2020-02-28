#ifndef _H_HOIBC_CONSTANTS
#define _H_HOIBC_CONSTANTS


#include "hoibc_types.hpp"
#include <cmath>

namespace hoibc
{
  const complex ci = complex(0.,1.); // imaginary unit
  const real    pi = 4.*std::atan(1.);
  const real    speed_of_light = 299792458; //  (m/s)
  const real    vacuum_impedance = 376.730452; // (Ohm)

  real free_space_impedance(const real& frequency);
}

#endif