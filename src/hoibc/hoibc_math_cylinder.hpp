#ifndef _H_HOIBC_MATH_CYLINDER
#define _H_HOIBC_MATH_CYLINDER

#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc {
  namespace cylinder {

    matrix<complex> JE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> HE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> JH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> HH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> MJ(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);
    matrix<complex> MH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);

    big_matrix<complex> impedance_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& inner_radius);

    big_matrix<complex> reflexion_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& inner_radius);
  }
}
#endif