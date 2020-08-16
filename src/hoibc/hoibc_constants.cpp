#include "hoibc_constants.hpp"

using namespace hoibc;

real hoibc::free_space_wavenumber(const real& frequency){
  return 2.*pi*frequency*1.E9/speed_of_light; // assuming the frequency is in GHz
};