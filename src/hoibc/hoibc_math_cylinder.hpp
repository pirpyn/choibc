#ifndef _H_HOIBC_MATH_CYLINDER
#define _H_HOIBC_MATH_CYLINDER

#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc {
  namespace cylinder {

    matrix<complex> AE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> BE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> AH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);
    matrix<complex> BH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar);

    matrix<complex> MA(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);
    matrix<complex> MB(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);
    matrix<complex> NE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);
    matrix<complex> NH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B);

    big_matrix<complex> impedance_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& inner_radius);
    big_matrix<complex> reflexion_from_impedance(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> impedance, const real& outer_radius);

    big_matrix<complex> reflexion_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius);
    big_matrix<complex> impedance_from_reflexion(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> reflexion, const real& outer_radius);

    void get_matrices_AB_EH(const real& radius, const array<real>& vn, const array<real>& vkz, const complex& k, const complex& etar, big_matrix<complex>& mJE, big_matrix<complex>& mHE, big_matrix<complex>& mJH, big_matrix<complex>& mHH);
    void get_matrices_LD_LR(const real& radius, const array<real>& vn, const array<real>& vkz, big_matrix<real>& LD = empty_bigmatrix_real, big_matrix<real>& LR = empty_bigmatrix_real);
  }
}
#endif