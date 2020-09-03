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

matrix<complex> hoibc::cylinder::NE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
  return JE(radius,n,kz,k,etar) + HE(radius,n,kz,k,etar)*B;
}

matrix<complex> hoibc::cylinder::NH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
  return JH(radius,n,kz,k,etar) + HH(radius,n,kz,k,etar)*B;
}

big_matrix<complex> hoibc::cylinder::impedance_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius){
  big_matrix<complex> impedance_ex = big_init(vn.size(),vkz.size(),material.initial_impedance);

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::cylinder::impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
    std::exit(1);
  }

  real radius = inner_radius;

  for (std::size_t l = 0; l < material.thickness.size() - 1; l++) {
    complex mu = material.mur[l];
    complex eps = material.epsr[l];
    complex etar, nur;
    check_and_set_material(l+1, eps, mu, etar, nur,material.loss);
    real d = material.thickness[l];

    const complex k = k0*nur;

    for (unsigned int i=0; i<vn.size();i++) {
      const real n = vn[i];
      for (unsigned int j=0; j<vkz.size();j++) {
        const real kz = vkz[j];

        const matrix<complex>& B = impedance_ex[i][j];

        const matrix<complex> mMJ = inv(MJ(radius,n,kz,k,etar,B));
        const matrix<complex> mMH = inv(MH(radius,n,kz,k,etar,B));

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

big_matrix<complex> hoibc::cylinder::reflexion_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius){
  big_matrix<complex> reflexion_ex = big_init(vn.size(),vkz.size(),complex(0.,0.));

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::cylinder::impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
    std::exit(1);
  }

  // The deepest interface ( i.e between pec and coating )
  complex mu_upper = material.mur[0];
  complex eps_upper = material.epsr[0];
  complex etar_upper, nur_upper;
  complex k_upper = k0*nur_upper;

  check_and_set_material(1, eps_upper, mu_upper, etar_upper, nur_upper, material.loss);

  real radius = inner_radius;

  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    for (unsigned int j=0; j<vkz.size();j++) {
      const real kz = vkz[j];

      reflexion_ex[i][j] = 
        - MH(radius,n,kz,k_upper,etar_upper,material.initial_impedance) 
        % MJ(radius,n,kz,k_upper,etar_upper,material.initial_impedance);
    }
  }

  // Strictly intermediate interfaces
  for (std::size_t l = 0; l < material.thickness.size()-1; l++) {
    complex mu_lower = mu_upper;
    complex eps_lower = mu_upper;
    complex etar_lower = etar_upper;
    complex nur_lower = nur_upper;

    complex k_lower = k_upper;

    mu_upper = material.mur[l+1];
    eps_upper = material.epsr[l+1];
    check_and_set_material(l+1, eps_upper, mu_upper, etar_upper, nur_upper, material.loss);

    k_upper = k0*nur_upper;

    radius += material.thickness[l];

    for (unsigned int i=0; i<vn.size();i++) {
      const real n = vn[i];
      for (unsigned int j=0; j<vkz.size();j++) {
        const real kz = vkz[j];

        const matrix<complex>& B = reflexion_ex[i][j];

        const matrix<complex> mNE = inv(NE(radius,n,kz,k_lower,etar_lower,B));
        const matrix<complex> mNH = inv(NH(radius,n,kz,k_lower,etar_lower,B));

        const matrix<complex> mJE = JE(radius,n,kz,k_upper,etar_upper);
        const matrix<complex> mHE = HE(radius,n,kz,k_upper,etar_upper);

        const matrix<complex> mJH = JH(radius,n,kz,k_upper,etar_upper);
        const matrix<complex> mHH = HH(radius,n,kz,k_upper,etar_upper);

        reflexion_ex[i][j] = - ( mNE*mHE - mNH*mHH ) % ( mNE*mJE - mNH*mJH );
      }
    }
  }

  // The coating-vacuum interface (radius = r_ext)
  // last layer
  complex mu_lower = mu_upper;
  complex eps_lower = mu_upper;
  complex etar_lower = etar_upper;
  complex nur_lower = nur_upper;
  complex k_lower = k_upper;

  // vacuum
  radius += material.thickness[material.thickness.size()-1];
  mu_upper = 1.;
  eps_upper = 1.;
  etar_upper = 1.;
  nur_upper = 1.;
  k_upper = k0;

  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    for (unsigned int j=0; j<vkz.size();j++) {
      const real kz = vkz[j];

      //if (std::abs(std::pow(k_upper,2) - std::pow(kz,2)) > 0.) {
        const matrix<complex>& B = reflexion_ex[i][j];

        const matrix<complex> mNE = inv(NE(radius,n,kz,k_lower,etar_lower,B));
        const matrix<complex> mNH = inv(NH(radius,n,kz,k_lower,etar_lower,B));

        const matrix<complex> mJE = JE(radius,n,kz,k_upper,etar_upper);
        const matrix<complex> mHE = HE(radius,n,kz,k_upper,etar_upper);

        const matrix<complex> mJH = JH(radius,n,kz,k_upper,etar_upper);
        const matrix<complex> mHH = HH(radius,n,kz,k_upper,etar_upper);

        reflexion_ex[i][j] = - ( mNE*mHE - mNH*mHH ) % ( mNE*mJE - mNH*mJH );

//      } else {
//        std::cerr << "warning: reflexion_infinite_cylinder: can't compute Hankel function at k3*outer_radius for k3 = 0." << std::endl;
//        R(n1,n2,:,:) = reshape(cmplx([-1.,0.,0.,1.],[0.,0.,0.,0.],wp),[2,2])
//      }
    }
  }
  return reflexion_ex;
}

big_matrix<complex> hoibc::cylinder::impedance_from_reflexion(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> reflexion, const real& outer_radius){
  big_matrix<complex> impedance = big_init(vn.size(),vkz.size(),complex(0.,0.));
  // In the vacuum
  complex k = k0;
  complex etar = 1.;
  real radius = outer_radius;
  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    for (unsigned int j=0; j<vkz.size();j++) {
      const real kz = vkz[j];
      if (std::abs(std::pow(k,2) - std::pow(kz,2)) > 0.) {
          matrix<complex> mR = reflexion[i][j];
          matrix<complex> mNE = NE(radius,n,kz,k,etar,mR);
          matrix<complex> mNH = NH(radius,n,kz,k,etar,mR);
          impedance[i][j] = mNE / mNH;
      } else {
          std::cerr << "warning: impedance_from_reflexion: can't compute Hankel function at k3*outer_radius for k3 = 0" << std::endl;
          // With big init we have impedance[i][j] = 0
      }
    }
  }
  return impedance;
}

void hoibc::cylinder::get_matrices_LD_LR(const real& radius, const array<real>& vn, const array<real>& vkz, big_matrix<real>& LD, big_matrix<real>& LR){
  LD.resize(0);
  LR.resize(0);
  LD = big_init(vn.size(),vkz.size(),0.);
  LR = big_init(vn.size(),vkz.size(),0.);

  for (std::size_t i = 0; i < vn.size(); i++){
    real n = vn[i];
    for (std::size_t j = 0; j < vkz.size(); j++){
      real kz = vkz[j];
      LD[i][j][0][0] = - std::pow(n/radius,2);
      LD[i][j][0][1] = - kz*n/radius;
      LD[i][j][1][0] =   LD[i][j][0][1];
      LD[i][j][1][1] = - std::pow(kz,2);

      LR[i][j][0][0] = - LD[i][j][1][1];
      LR[i][j][0][1] =   LD[i][j][0][1];
      LR[i][j][1][0] =   LD[i][j][0][1];
      LR[i][j][1][1] = - LD[i][j][0][0];
    }
  }
}

void hoibc::cylinder::get_matrices_JH_EH(const real& radius, const array<real>& vn, const array<real>& vkz, const complex& k, const complex& etar, big_matrix<complex>& mJE ,big_matrix<complex>& mHE ,big_matrix<complex>& mJH ,big_matrix<complex>& mHH){
  for (std::size_t i = 0; i < vn.size(); i++){
    real n = vn[i];
    for (std::size_t j = 0; j < vkz.size(); j++){
      real kz = vkz[j];
      mJE[i][j] = JE(radius,vn[i],vkz[j],k,etar);
      mJH[i][j] = JH(radius,vn[i],vkz[j],k,etar);
      mHE[i][j] = HE(radius,vn[i],vkz[j],k,etar);
      mHH[i][j] = HH(radius,vn[i],vkz[j],k,etar);
    }
  }
}
