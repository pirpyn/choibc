#ifndef _H_HOIBC_CONSTANTS
#define _H_HOIBC_CONSTANTS

#include "hoibc_types.hpp"
#include <cmath>

namespace hoibc
{
  const complex ci = complex(0.,1.); // imaginary unit
  const real    pi = 3.141592653589793;
  const real    speed_of_light = 299792458; //  (m/s)
  const real    vacuum_impedance = 376.730452; // (Ohm)
  const real    sqrt_two = 1.4142135623730951;

  real free_space_wavenumber(const real& frequency);
}

#endif