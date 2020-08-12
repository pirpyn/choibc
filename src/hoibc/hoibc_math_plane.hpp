#ifndef _H_HOIBC_MATH_PLANE
#define _H_HOIBC_MATH_PLANE

#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc {
  namespace plane {
    matrix<complex> AE(const real& kx, const real& ky, const complex& k, const real& z);
    matrix<complex> BE(const real& kx, const real& ky, const complex& k, const real& z);
    matrix<complex> AH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z);
    matrix<complex> BH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z);
    void get_matrices_AB_EH(const std::vector<real>& f1, const std::vector<real>& f2, const complex& k, const complex& etar, const real& z, big_matrix<complex>& mAE ,big_matrix<complex>& mBE ,big_matrix<complex>& mAH ,big_matrix<complex>& mBH);
  }
  matrix<complex> MA(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp);
  matrix<complex> MB(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp);

  big_matrix<complex> impedance_infinite_plane(const std::vector<real>& vkx, const std::vector<real>& vky, const real& k0, const material_t& material);
  big_matrix<complex> reflexion_infinite_plane(const std::vector<real>& vkx, const std::vector<real>& vky, const real& k0, const material_t& material);
}
#endif