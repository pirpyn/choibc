#include "hoibc_math_cylinder.hpp"
#include "hoibc_math.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_bessel.hpp"
#include <cmath> // std:: pow, std::sqrt
#include <iostream>

using namespace hoibc;

  // The tangential part of E / n x H writes in each layer where C1,C2 depends on the layer
  // E_t(r,theta,z)     = 1/(2 pi) sum_n int_R JE(r,n,kz) C1(kz,n) + HE(r,n,kz) C2(kz,n) dkz
  // n x H_t(r,theta,z) = 1/(2 pi) sum_n int_R JH(r,n,kz) C1(kz,n) + HH(r,n,kz) C2(kz,n) dkz

  // See hoibc_math_plane.cpp for more

matrix<complex> hoibc::cylinder::JE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
  matrix<complex> JE;
  const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
  JE[0][0] = -n*kz/(radius*std::pow(k3,2))*bessel1(n,k3*radius);
  JE[1][0] = bessel1(n,k3*radius);
  JE[0][1] = ci*etar*k/k3*bessel1p(n,k3*radius);
  JE[1][1] = complex(0.,0.);
  return JE;
}

matrix<complex> hoibc::cylinder::HE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
  matrix<complex> HE;
  const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
  HE[0][0] = -n*kz/(radius*std::pow(k3,2))*bessel2(n,k3*radius);
  HE[1][0] = bessel2(n,k3*radius);
  HE[0][1] = ci*etar*k/k3*bessel2p(n,k3*radius);
  HE[1][1] = complex(0.,0.);
  return HE;
}

matrix<complex> hoibc::cylinder::JH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
  matrix<complex> JH;
  const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
  JH[0][0] = complex(0.,0.);
  JH[1][0] = -ci*etar*k/k3*bessel1p(n,k3*radius);
  JH[0][1] = -bessel1(n,k3*radius);
  JH[1][1] = -n*kz/(radius*std::pow(k3,2))*bessel1(n,k3*radius);
  return JH;
}

matrix<complex> hoibc::cylinder::HH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
  matrix<complex> HH;
  const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
  HH[0][0] = complex(0.,0.);
  HH[1][0] = -ci*etar*k/k3*bessel2p(n,k3*radius);
  HH[0][1] = -bessel2(n,k3*radius);
  HH[1][1] = -n*kz/(radius*std::pow(k3,2))*bessel2(n,k3*radius);
  return HH;
}

matrix<complex> hoibc::cylinder::MJ(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
  return JE(radius,n,kz,k,etar) - B*JH(radius,n,kz,k,etar);
}

matrix<complex> hoibc::cylinder::MH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
  return HE(radius,n,kz,k,etar) - B*HH(radius,n,kz,k,etar);
}

big_matrix<complex> hoibc::cylinder::impedance_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius){
  big_matrix<complex> impedance_ex = big_init(vn.size(),vkz.size(),material.initial_impedance);

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::cylinder::impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
    std::exit(1);
  }

  real radius = inner_radius;

  for (std::size_t i = 0; i < material.thickness.size(); i++) {
    complex mu = material.mur[i];
    complex eps = material.epsr[i];
    complex etar, nur;
    check_and_set_material(i+1, eps, mu, etar, nur,material.loss);
    real d = material.thickness[i];

    const complex k = k0*nur;

    for (unsigned int i=0; i<vn.size();i++) {
      const real n = vn[i];
      for (unsigned int j=0; j<vkz.size();j++) {
        const real kz = vkz[j];

        const matrix<complex>& B = impedance_ex[i][j];

        const matrix<complex> mMJ = inv<complex>(MJ(radius,n,kz,k,etar,B));
        const matrix<complex> mMH = inv<complex>(MH(radius,n,kz,k,etar,B));

        const matrix<complex> mJE = JE(radius+d,n,kz,k,etar);
        const matrix<complex> mHE = HE(radius+d,n,kz,k,etar);

        const matrix<complex> mJH = JH(radius+d,n,kz,k,etar);
        const matrix<complex> mHH = HH(radius+d,n,kz,k,etar);

        // Using the overloaded / operator A / B = A*B^{-1}
        impedance_ex[i][j] = (mHE*mMH - mJE*mMJ) / (mHH*mMH - mJH*mMJ);

      }
    }
    radius = radius + d;
  }
  return impedance_ex;
}


big_matrix<complex> hoibc::cylinder::reflexion_from_impedance(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> impedance, const real& outer_radius){
  big_matrix<complex> reflexion = big_init(vn.size(),vkz.size(),complex(0.,0.));
  // in the vacuum
  real radius = outer_radius;
  complex k = k0;
  complex etar = 1.;
  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    for (unsigned int j=0; j<vkz.size();j++) {
      const real kz = vkz[j];
      if (std::abs(std::pow(k,2) - std::pow(kz,2)) <= std::numeric_limits<real>::epsilon()){
        // In that case, we get a 0/0 which result in NaN.
        // Mathematical analysis gives to following value at that point
        matrix<complex> mR0 { 0. };
        mR0[0][0] = complex(-1.,0.);
        mR0[1][1] = complex(1.,0.);
        reflexion[i][j] = mR0;
      } else {
        //Using the overloaded % operator A % B = A^{-1}*B
        reflexion[i][j] = - MH(radius,n,kz,k,etar,impedance[i][j]) % MJ(radius,n,kz,k,etar,impedance[i][j]);
      }
    }
  }
  return reflexion;
}

big_matrix<complex> hoibc::cylinder::reflexion_infinite(const array<real>& n, const array<real>& kz, const real& k0, const material_t& material, const real& radius){
  big_matrix<complex> reflexion_ex;
  return reflexion_ex;
  NOTFINISHED("hoibc_math_cylinder::reflexion_infinite");
}
