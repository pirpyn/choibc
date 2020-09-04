#include "hoibc_math_sphere.hpp"
#include "hoibc_math.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_bessel.hpp"
#include <cmath> // std:: pow, std::sqrt
#include <iostream>

using namespace hoibc;

// The tangential part of E or n x H writes in each layer, where C1/C2 depends on the layer
// E_t(r,theta,phi)     =  sum_n sum_m D(n,m,theta,phi) ( AE(r,n,m) C1(n,m) + BE(r,n,m) C2(n,m) )
// n x H_t(r,theta,phi) =  sum_n sum_m D(n,m,theta,phi) ( AH(r,n,m) C1(n,m) + BH(r,n,m) C2(n,m) )

// The matrix D includes the spherical harmonics functions Pnm(cos(theta))exp(i m phi)
// but is not present in the impedance

// See mod_hoibc_math_plane.cpp for more

// the 'tilde' bessel functions, i.e d/dx (xz_n(x)) where z_n is a spherical bessel or hankel
static complex sbessel1t(const real& nu, const complex& z){
  return sbessel1(nu,z) + z*sbessel1p(nu,z);
}

static complex sbessel2t(const real& nu, const complex& z){
  return sbessel2(nu,z) + z*sbessel2p(nu,z);
}

matrix<complex> hoibc::sphere::AE(const real& radius,const real& n, const complex& k){
  matrix<complex> AE;
  AE[0][0] = sbessel1(n,k*radius);
  AE[1][0] = complex(0.,0.);
  AE[0][1] = complex(0.,0.);
  AE[1][1] = sbessel1t(n,k*radius);
  return AE;
}

matrix<complex> hoibc::sphere::BE(const real& radius,const real& n, const complex& k){
  matrix<complex> BE;
  BE[0][0] = sbessel2(n,k*radius);
  BE[1][0] = complex(0.,0.);
  BE[0][1] = complex(0.,0.);
  BE[1][1] = sbessel2t(n,k*radius);
  return BE;
}

matrix<complex> hoibc::sphere::AH(const real& radius,const real& n, const complex& k, const complex& etar){
  matrix<complex> AH;
  AH[0][0] = -ci/etar/k/radius*sbessel1t(n,k*radius);
  AH[1][0] = complex(0.,0.);
  AH[0][1] = complex(0.,0.);
  AH[1][1] = ci/etar*k*radius*sbessel1(n,k*radius);
  return AH;
}

matrix<complex> hoibc::sphere::BH(const real& radius,const real& n, const complex& k, const complex& etar){
  matrix<complex> BH;
  BH[0][0] = -ci/etar/k/radius*sbessel2t(n,k*radius);
  BH[1][0] = complex(0.,0.);
  BH[0][1] = complex(0.,0.);
  BH[1][1] = ci/etar*k*radius*sbessel2(n,k*radius);
  return BH;
}

matrix<complex> hoibc::sphere::MA(const real& radius, const real& n, const complex& k, const complex& etar, const matrix<complex>& B){
  return AE(radius,n,k) - B*AH(radius,n,k,etar);
}

matrix<complex> hoibc::sphere::MB(const real& radius, const real& n, const complex& k, const complex& etar, const matrix<complex>& B){
  return BE(radius,n,k) - B*BH(radius,n,k,etar);
}

matrix<complex> hoibc::sphere::NE(const real& radius, const real& n, const complex& k, const matrix<complex>& B){
  return AE(radius,n,k) + BE(radius,n,k)*B;
}

matrix<complex> hoibc::sphere::NH(const real& radius, const real& n, const complex& k, const complex& etar, const matrix<complex>& B){
  return AH(radius,n,k,etar) + BH(radius,n,k,etar)*B;
}

big_matrix<complex> hoibc::sphere::impedance_infinite(const array<real>& vn, const real& k0, const material_t& material, const real& inner_radius){
  big_matrix<complex> impedance_ex = big_init(vn.size(),1,material.initial_impedance);

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::sphere::impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
    std::exit(1);
  }

  real radius = inner_radius;

  for (std::size_t l = 0; l < material.thickness.size(); l++) {
    complex mu = material.mur[l];
    complex eps = material.epsr[l];
    complex etar, nur;
    check_and_set_material(l+1, eps, mu, etar, nur,material.loss);
    real d = material.thickness[l];

    const complex k = k0*nur;

    for (unsigned int i=0; i<vn.size();i++) {
      const real n = vn[i];

        const matrix<complex>& B = impedance_ex[i][0];

        const matrix<complex> mMA = inv(MA(radius,n,k,etar,B));
        const matrix<complex> mMB = inv(MB(radius,n,k,etar,B));

        const matrix<complex> mAE = AE(radius+d,n,k);
        const matrix<complex> mBE = BE(radius+d,n,k);

        const matrix<complex> mAH = AH(radius+d,n,k,etar);
        const matrix<complex> mBH = BH(radius+d,n,k,etar);

        // Using the overloaded / operator A / B = A*B^{-1}
        impedance_ex[i][0] = (mBE*mMB - mAE*mMA) / (mBH*mMB - mAH*mMA);
      }
    radius = radius + d;
  }
  return impedance_ex;
}


big_matrix<complex> hoibc::sphere::reflexion_from_impedance(const array<real>& vn, const real& k0, const big_matrix<complex> impedance, const real& outer_radius){
  big_matrix<complex> reflexion = big_init(vn.size(),1,complex(0.,0.));
  // in the vacuum
  real radius = outer_radius;
  complex k = k0;
  complex etar = 1.;
  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    const matrix<complex>& B = impedance[i][1];
    reflexion[i][0] = - MB(radius,n,k,etar,B) % MA(radius,n,k,etar,B);
  }
  return reflexion;
}

big_matrix<complex> hoibc::sphere::reflexion_infinite(const array<real>& vn, const real& k0, const material_t& material, const real& inner_radius){
  big_matrix<complex> reflexion_ex = big_init(vn.size(),1,complex(0.,0.));

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::sphere::impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
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

    reflexion_ex[i][1] = 
      - MB(radius,n,k_upper,etar_upper,material.initial_impedance) 
      % MA(radius,n,k_upper,etar_upper,material.initial_impedance);
  }

  // Strictly intermediate interfaces
  for (std::size_t l = 0; l < material.thickness.size()-1; l++) {
    complex etar_lower = etar_upper;
    complex k_lower = k_upper;

    mu_upper = material.mur[l+1];
    eps_upper = material.epsr[l+1];
    check_and_set_material(l+1, eps_upper, mu_upper, etar_upper, nur_upper, material.loss);

    k_upper = k0*nur_upper;

    radius += material.thickness[l];

    for (unsigned int i=0; i<vn.size();i++) {
      const real n = vn[i];

      const matrix<complex>& B = reflexion_ex[i][1];

      const matrix<complex> mNE = inv(NE(radius,n,k_lower,B));
      const matrix<complex> mNH = inv(NH(radius,n,k_lower,etar_lower,B));

      const matrix<complex> mAE = AE(radius,n,k_upper);
      const matrix<complex> mBE = BE(radius,n,k_upper);

      const matrix<complex> mAH = AH(radius,n,k_upper,etar_upper);
      const matrix<complex> mBH = BH(radius,n,k_upper,etar_upper);

      reflexion_ex[i][1] = - ( mNE*mBE - mNH*mBH ) % ( mNE*mAE - mNH*mAH );
    }
  }

  // The coating-vacuum interface (radius = r_ext)
  // last layer
  complex etar_lower = etar_upper;
  complex k_lower = k_upper;

  // vacuum
  radius += material.thickness[material.thickness.size()-1];
  etar_upper = 1.;
  k_upper = k0;

  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];

    const matrix<complex>& B = reflexion_ex[i][1];

    const matrix<complex> mNE = inv(NE(radius,n,k_lower,B));
    const matrix<complex> mNH = inv(NH(radius,n,k_lower,etar_lower,B));

    const matrix<complex> mAE = AE(radius,n,k_upper);
    const matrix<complex> mBE = BE(radius,n,k_upper);

    const matrix<complex> mAH = AH(radius,n,k_upper,etar_upper);
    const matrix<complex> mBH = BH(radius,n,k_upper,etar_upper);

    reflexion_ex[i][1] = - ( mNE*mBE - mNH*mBH ) % ( mNE*mAE - mNH*mAH );

  }
  return reflexion_ex;
}

big_matrix<complex> hoibc::sphere::impedance_from_reflexion(const array<real>& vn, const real& k0, const big_matrix<complex> reflexion, const real& outer_radius){
  big_matrix<complex> impedance = big_init(vn.size(),1,complex(0.,0.));
  // In the vacuum
  complex k = k0;
  complex etar = 1.;
  real radius = outer_radius;
  for (unsigned int i=0; i<vn.size();i++) {
    const real n = vn[i];
    const matrix<complex> mR = reflexion[i][1];
    const matrix<complex> mNE = NE(radius,n,k,mR);
    const matrix<complex> mNH = NH(radius,n,k,etar,mR);
    impedance[i][1] = mNE / mNH;
    }
  return impedance;
}

void hoibc::sphere::get_matrices_LD_LR(const real& radius, const array<real>& vn, big_matrix<real>& LD, big_matrix<real>& LR){
  LD.resize(0);
  LR.resize(0);
  LD = big_init(vn.size(),1,0.);
  LR = big_init(vn.size(),1,0.);

  for (std::size_t i = 0; i < vn.size(); i++){
    real n = vn[i];
    LD[i][1][1][1] = - n*(n+1)/std::pow(radius,2);

    LR[i][1][0][0] = n*(n+1)/std::pow(radius,2);

  }
}

void hoibc::sphere::get_matrices_AB_EH(const real& radius, const array<real>& vn, const complex& k, const complex& etar, big_matrix<complex>& mAE ,big_matrix<complex>& mBE ,big_matrix<complex>& mAH ,big_matrix<complex>& mBH){
  mAE = big_init(vn.size(),1,complex(0.,0.));
  mBE = big_init(vn.size(),1,complex(0.,0.));
  mAH = big_init(vn.size(),1,complex(0.,0.));
  mBH = big_init(vn.size(),1,complex(0.,0.));
  for (std::size_t i = 0; i < vn.size(); i++){
    real n = vn[i];
    mAE[i][1] = AE(radius,n,k);
    mAH[i][1] = AH(radius,n,k,etar);
    mBE[i][1] = BE(radius,n,k);
    mBH[i][1] = BH(radius,n,k,etar);
  }
}
