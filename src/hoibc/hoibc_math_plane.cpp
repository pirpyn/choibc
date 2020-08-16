#include "hoibc_math_plane.hpp"
#include "hoibc_types.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_math.hpp" // matmul, inv, big_init, overloaded *+-/
#include <iostream>
#include <numeric> // std::accumulate
#include <cmath> // pow, sqrt, imag, real, abs
#include <limits> // epsilon

using namespace hoibc;
using namespace hoibc::plane;

matrix<complex> hoibc::plane::AE(const real& kx, const real& ky, const complex& k, const real& z){
  complex k3;
  if ((z>=0.)&&(std::abs(std::imag(k))<=0.)) {
    // On the outside  where z >= 0 and k = k0 is real, we need non-infinite solution 
    // when z goes to infinity hence we need the imaginary part of k3 to be negative 
    // for exp(-i k3 z) to go to 0 when z goes to infinity.
    // That's why there is a minus in front of the zero below
    k3 = std::sqrt(complex(std::pow(std::real(k),2) - std::pow(kx,2) - std::pow(ky,2),-0.));
  }
  else {
    k3 = std::sqrt(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2));
  }
  matrix<complex> AE;
  AE[0][0] = std::exp(ci*k3*z)*ci*k3;
  AE[1][0] = 0.;
  AE[0][1] = AE[1][0];
  AE[1][1] = AE[0][0];
  return AE;
}

matrix<complex> hoibc::plane::BE(const real& kx, const real& ky, const complex& k, const real& z){
  complex k3;
  if ((z>=0.)&&(std::abs(std::imag(k))<=0.)) {
    // See AE above
    k3 = std::sqrt(complex(std::pow(std::real(k),2) - std::pow(kx,2) - std::pow(ky,2),-0.));
  }
  else {
    k3 = std::sqrt(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2));
  }
  matrix<complex> BE;
  BE[0][0] = -std::exp(-ci*k3*z)*ci*k3;
  BE[1][0] = 0.;
  BE[0][1] = BE[1][0];
  BE[1][1] = BE[0][0];
  return BE;
}

matrix<complex> hoibc::plane::AH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z){
  complex k3;
  if ((z>=0.)&&(std::abs(std::imag(k))<=0.)) {
    // See AE above
    k3 = std::sqrt(complex(std::pow(std::real(k),2) - std::pow(kx,2) - std::pow(ky,2),-0.));
  }
  else {
    k3 = std::sqrt(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2));
  }
  matrix<complex> AH;
  AH[0][0] = std::exp(ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(ky,2));
  AH[1][0] = std::exp(ci*k3*z)*ci/k/etar*kx*ky;
  AH[0][1] = AH[1][0];
  AH[1][1] = std::exp(ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(kx,2));
  return AH;
}

matrix<complex> hoibc::plane::BH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z){
  complex k3;
  if ((z>=0.)&&(std::abs(std::imag(k))<=0.)) {
    // See AE above
    k3 = std::sqrt(complex(std::pow(std::real(k),2) - std::pow(kx,2) - std::pow(ky,2),-0.));
  }
  else {
    k3 = std::sqrt(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2));
  }
  matrix<complex> BH;
  BH[0][0] = std::exp(-ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(ky,2));
  BH[1][0] = std::exp(-ci*k3*z)*ci/k/etar*kx*ky;
  BH[0][1] = BH[1][0];
  BH[1][1] = std::exp(-ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(kx,2));
  return BH;
}

void hoibc::plane::get_matrices_AB_EH(const std::vector<real>& f1, const std::vector<real>& f2, const complex& k, const complex& etar, const real& z, big_matrix<complex>& mAE ,big_matrix<complex>& mBE ,big_matrix<complex>& mAH ,big_matrix<complex>& mBH){
  // return big 4 rank arrays that store AE, BE, AH, BH matrices
  // for every kx and ky
  //
  // e.g.
  //
  //    mAE[i][j] <=> AE(kx(i),ky(j),k,z)

  mAE = big_init(f1.size(),f2.size(),complex(0.,0.));
  mBE = big_init(f1.size(),f2.size(),complex(0.,0.));
  mAH = big_init(f1.size(),f2.size(),complex(0.,0.));
  mBH = big_init(f1.size(),f2.size(),complex(0.,0.));

  for (std::size_t i1 = 0; i1 < f1.size(); i1++){
    for (std::size_t i2 = 0; i2 < f2.size(); i2++){
      mAE[i1][i2] = AE(f1[i1],f2[i2],k,z);
      mBE[i1][i2] = BE(f1[i1],f2[i2],k,z);
      mAH[i1][i2] = AH(f1[i1],f2[i2],k,etar,z);
      mBH[i1][i2] = BH(f1[i1],f2[i2],k,etar,z);
    }
  }
}


matrix<complex> hoibc::MA(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp){
  return AE(kx,ky,k,z) - matmul(imp,AH(kx,ky,k,etar,z));
}

matrix<complex> hoibc::MB(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp){
  return BE(kx,ky,k,z) - matmul(imp,BH(kx,ky,k,etar,z));
}

big_matrix<complex> hoibc::plane::impedance_infinite(const std::vector<real> &vkx, const  std::vector<real> &vky, const real& k0, const material_t& material){

  big_matrix<complex> impedance_ex = big_init(vkx.size(), vky.size(), material.initial_impedance);

  const std::vector<real>& thickness = material.thickness;

  real h = - std::accumulate(thickness.begin(),thickness.end(),0.);

  for (std::size_t i = 0; i < thickness.size(); i++) {
    complex mu = material.mur[i];
    complex eps = material.epsr[i];
    real d = thickness[i];

    if ( (std::imag(mu)==0.) && (std::imag(eps)==0.) ) {
      std::cerr << "Layer " << i+1 << ": adding artificial loss of" << material.loss << std::endl;
      mu = complex(std::real(mu),-std::abs(material.loss));
      eps = complex(std::real(eps),-std::abs(material.loss));
    }

    const complex etar = std::sqrt(mu/eps);
    const complex nur = std::sqrt(mu*eps);

    if (std::imag(nur)>0.) {
      std::cerr << "error: hoibc::impedance_infinite_plane: Im(nur) > 0" << std::endl;
      exit(1);
    }
    if (std::real(etar)<0.) {
      std::cerr << "error: hoibc::impedance_infinite_plane: Re(etar) < 0" << std::endl;
      exit(1);
    }

    const complex k = k0*nur;

    for (unsigned int i=0; i<vkx.size();i++) {
      const real kx = vkx[i];
      for (unsigned int j=0; j<vky.size();j++) {
        const real ky = vky[j];

        const matrix<complex>& B = impedance_ex[i][j];

        const matrix<complex> mMB = inv<complex>(MB(kx,ky,k,etar,h,B));
        const matrix<complex> mMA = inv<complex>(MA(kx,ky,k,etar,h,B));

        const matrix<complex> mAE = AE(kx,ky,k,h+d);
        const matrix<complex> mBE = BE(kx,ky,k,h+d);

        const matrix<complex> mAH = AH(kx,ky,k,etar,h+d);
        const matrix<complex> mBH = BH(kx,ky,k,etar,h+d);

        // Using the overloaded / operator A / B = A*B^{-1}
        impedance_ex[i][j] = (mBE*mMB - mAE*mMA) / (mBH*mMB - mAH*mMA);
      }
    }
    h+=d;
  }
  return impedance_ex;
}

big_matrix<complex> hoibc::plane::impedance_from_reflexion(const std::vector<real>& vkx,const std::vector<real>& vky, const real& k0,const big_matrix<complex>& reflexion){
  // In the vacuum
  const complex k = k0;
  const complex etar = 1.;
  const real h = 0.;
  
  big_matrix<complex> impedance = big_init(vkx.size(),vky.size(),complex(0.,0.));

  for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
    const real kx = vkx[n1];
    for (std::size_t n2 = 0; n2 < vky.size(); n2++){
      const real ky = vky[n2];

      const matrix<complex> mAE = AE(kx,ky,k,h);
      const matrix<complex> mBE = BE(kx,ky,k,h);

      const matrix<complex> mAH = AH(kx,ky,k,etar,h);
      const matrix<complex> mBH = BH(kx,ky,k,etar,h);
      
      // Using the overloaded % operator A % B = A^{-1}*B
      impedance[n1][n2] = (mAH + mBH*reflexion[n1][n2]) % (mAE + mBE*reflexion[n1][n2]);
    }
  }
  return impedance;
}

big_matrix<complex> hoibc::plane::reflexion_infinite(const std::vector<real>& vkx, const std::vector<real>& vky, const real& k0, const material_t& material){

  big_matrix<complex> reflexion_ex = big_init(vkx.size(),vky.size(),complex(0.,0.));

  if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
    std::cerr << "error: hoibc::plane::reflexion_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
    std::exit(1);
  }

  // The deepest interface ( i.e between pec and coating )
  complex mu_upper = material.mur[0];
  complex eps_upper = material.epsr[0];
  if ((std::abs(std::imag(mu_upper))<=0.)&&(std::abs(std::imag(eps_upper))<=0.)&&(std::abs(material.loss)>0.)){
    std::cerr << "Layer " << 1 << ": adding artificial loss of" << material.loss << std::endl;
    mu_upper  = complex(std::real(mu_upper),-std::abs(material.loss));
    eps_upper = complex(std::real(eps_upper),-std::abs(material.loss));
  }

  complex etar_upper = std::sqrt(mu_upper/eps_upper);
  complex nur_upper = std::sqrt(mu_upper*eps_upper);

  if (std::imag(nur_upper) > 0.){
    std::cerr << "error: hoibc::plane::reflexion_infinite: Im(nur) > 0 (" << std::imag(nur_upper) << ")" << std::endl;
    exit(1);
  }
  if (std::real(etar_upper) < 0.){
    std::cerr <<"error: hoibc::plane::reflexion_infinite: Re(etar) < 0 (" << std::real(etar_upper) << ")" << std::endl;
    exit(1);
  }

  complex k_upper = k0*nur_upper;

  real h = - std::accumulate(material.thickness.begin(),material.thickness.end(),0.);

  for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
    const real kx = vkx[n1];
    for (std::size_t n2 = 0; n2 < vky.size(); n2++){
      const real ky = vky[n2];

      const matrix<complex> mAE = AE(kx,ky,k_upper,h);
      const matrix<complex> mBE = BE(kx,ky,k_upper,h);

      const matrix<complex> mAH = AH(kx,ky,k_upper,etar_upper,h);
      const matrix<complex> mBH = BH(kx,ky,k_upper,etar_upper,h);

      reflexion_ex[n1][n2] = - (mBE - material.initial_impedance*mBH)%(mAE - material.initial_impedance*mAH);
    }
  }

  // Strictly intermediate interfaces ( i.e label 2 to n-1 )
  for (std::size_t l = 0; l < material.thickness.size() - 1; l++){
    complex mu_lower = material.mur[l];
    complex eps_lower = material.epsr[l];
    if ((std::abs(std::imag(mu_lower))<=0.)&&(std::abs(std::imag(eps_lower))<=0.)){
      mu_lower  = complex(std::real(mu_lower),-std::abs(material.loss));
      eps_lower = complex(std::real(eps_lower),-std::abs(material.loss));
    }
    complex etar_lower = std::sqrt(mu_lower/eps_lower);
    complex nur_lower = std::sqrt(mu_lower*eps_lower);
    complex k_lower = k0*nur_lower;

    mu_upper = material.mur[l+1];
    eps_upper = material.epsr[l+1];
    if ((std::abs(std::imag(mu_upper))<=0.)&&(std::abs(std::imag(eps_upper))<=0.)&&(std::abs(material.loss)>0.)){
      std::cerr << "Layer " << l+1 << ": adding artificial loss of" << material.loss << std::endl;
      mu_upper  = complex(std::real(mu_upper),-std::abs(material.loss));
      eps_upper = complex(std::real(eps_upper),-std::abs(material.loss));
    }
    etar_upper = std::sqrt(mu_upper/eps_upper);
    nur_upper = std::sqrt(mu_upper*eps_upper);
    k_upper = k0*nur_upper;

    if (std::imag(nur_upper) > 0.){
      std::cerr << "error: reflexion_infinite_plane: Im(nur) > 0" << std::endl;
      exit(1);
    }
    if (std::real(etar_upper) < 0.){
      std::cerr << "error: impedance_infinite_plane: Re(etar) < 0" << std::endl;
    }

    h += material.thickness[l];
    for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
      real kx = vkx[n1];
      for (std::size_t n2 = 0; n2 < vky.size(); n2++){
        real ky = vky[n2];
        
        matrix<complex>& mR = reflexion_ex[n1][n2];
        const matrix<complex> mAE_lower = AE(kx,ky,k_lower,h);
        const matrix<complex> mBE_lower = BE(kx,ky,k_lower,h);
        const matrix<complex> mAE_upper = AE(kx,ky,k_upper,h);
        const matrix<complex> mBE_upper = BE(kx,ky,k_upper,h);

        const matrix<complex> mAH_lower = AH(kx,ky,k_lower,etar_lower,h);
        const matrix<complex> mBH_lower = BH(kx,ky,k_lower,etar_lower,h);
        const matrix<complex> mAH_upper = AH(kx,ky,k_upper,etar_upper,h);
        const matrix<complex> mBH_upper = BH(kx,ky,k_upper,etar_upper,h);

        // Using the custom % operator A % B = A^(-1)*B
        reflexion_ex[n1][n2] = 
          - ((mAE_lower + mBE_lower*mR)%mBE_upper - (mAH_lower + mBH_lower*mR)%mBH_upper)
          % ((mAE_lower + mBE_lower*mR)%mAE_upper - (mAH_lower + mBH_lower*mR)%mAH_upper);
      }
    }
  }
  // The coating-vacuum interface (h = 0)
  // last layer
  complex mu_lower = material.mur[material.mur.size()-1];
  complex eps_lower = material.epsr[material.epsr.size()-1];
  if ((std::abs(std::imag(mu_lower))<=0.)&&(std::abs(std::imag(eps_lower))<=0.)){
    mu_lower  = complex(std::real(mu_lower),-std::abs(material.loss));
    eps_lower = complex(std::real(eps_lower),-std::abs(material.loss));
  }
  complex etar_lower = std::sqrt(mu_lower/eps_lower);
  complex nur_lower = std::sqrt(mu_lower*eps_lower);
  complex k_lower = k0*nur_lower;

  // vacuum
  h = 0.;
  mu_upper = 1.;
  eps_upper = 1.;
  etar_upper = 1.;
  nur_upper = 1.;
  k_upper = k0;

  for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
    real kx = vkx[n1];
    for (std::size_t n2 = 0; n2 < vky.size(); n2++){
      real ky = vky[n2];
      if (std::abs(std::pow(k_upper,2) - std::pow(kx,2) - std::pow(ky,2)) <= std::numeric_limits<real>::epsilon()){
        // In that case, we get a 0/0 which result in NaN.
        // Mathematical analysis gives to following value at that point
        matrix<complex> mR0 { 0. };
        mR0[0][0] = complex(1.,0.);
        mR0[1][1] = complex(-1.,0.);
        reflexion_ex[n1][n2] = mR0;
      } else {
        matrix<complex>& mR = reflexion_ex[n1][n2];
        const matrix<complex> mAE_lower = AE(kx,ky,k_lower,h);
        const matrix<complex> mBE_lower = BE(kx,ky,k_lower,h);
        const matrix<complex> mAE_upper = AE(kx,ky,k_upper,h);
        const matrix<complex> mBE_upper = BE(kx,ky,k_upper,h);

        const matrix<complex> mAH_lower = AH(kx,ky,k_lower,etar_lower,h);
        const matrix<complex> mBH_lower = BH(kx,ky,k_lower,etar_lower,h);
        const matrix<complex> mAH_upper = AH(kx,ky,k_upper,etar_upper,h);
        const matrix<complex> mBH_upper = BH(kx,ky,k_upper,etar_upper,h);

        reflexion_ex[n1][n2] =
            - ((mAE_lower + mBE_lower*mR)%mBE_upper - (mAH_lower + mBH_lower*mR)%mBH_upper)
            % ((mAE_lower + mBE_lower*mR)%mAE_upper - (mAH_lower + mBH_lower*mR)%mAH_upper);
      }
    }
  }
  return reflexion_ex;
}

big_matrix<complex> hoibc::plane::reflexion_from_impedance(const std::vector<real>& vkx,const std::vector<real>& vky,const real& k0, const big_matrix<complex>& impedance){
  
  big_matrix<complex> reflexion = big_init(vkx.size(),vky.size(),complex(0.,0.));

  // in the vacuum
  const real h = 0.;
  const complex k = k0;
  const complex etar = 1.;
  for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
    real kx = vkx[n1];
    for (std::size_t n2 = 0; n2 < vky.size(); n2++){
      real ky = vky[n2];
      if (std::abs(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2)) <= std::numeric_limits<real>::epsilon()){
        // In that case, we get a 0/0 which result in NaN.
        // Mathematical analysis gives to following value at that point
        matrix<complex> mR0 { 0. };
        mR0[0][0] = complex(1.,0.);
        mR0[1][1] = complex(1.,0.);
        reflexion[n1][n2] = mR0;
        if (std::abs(kx)<=0.){
          mR0[0][0] = complex(1.,0.);
          mR0[1][1] = complex(-1.,0.);
          reflexion[n1][n2] = mR0;
        }
        if (std::abs(ky)<=0.){
          mR0[0][0] = complex(-1.,0.);
          mR0[1][1] = complex(1.,0.);
          reflexion[n1][n2] = mR0;
        }
      } else {
        const matrix<complex> mZ  = impedance[n1][n2];
        const matrix<complex> mAE = AE(kx,ky,k,h);
        const matrix<complex> mBE = BE(kx,ky,k,h);
        const matrix<complex> mAH = AH(kx,ky,k,etar,h);
        const matrix<complex> mBH = BH(kx,ky,k,etar,h);

        reflexion[n1][n2] = - (mBE - mZ*mBH)%(mAE - mZ*mAH);
      }
    }
  }
  return reflexion;
}

void hoibc::plane::get_matrices_LD_LR(const std::vector<real>& vkx, const std::vector<real>& vky, big_matrix<real>& LD, big_matrix<real>& LR){
  LD.clear();
  LR.clear();
  LD = big_init(vkx.size(),vky.size(),0.);
  LR = big_init(vkx.size(),vky.size(),0.);

  // TODO Possible improve perf by doing
  // for n1
  //    for n2 
  //        tmp_vect.push_pack(tmp_mat)
  //    A.push_back(tmp_vect)

  for (std::size_t n1 = 0; n1 < vkx.size(); n1++){
    real kx = vkx[n1];
    for (std::size_t n2 = 0; n2 < vky.size(); n2++){
      real ky = vky[n2];
      LD[n1][n2][0][0] = -kx*kx;
      LD[n1][n2][0][1] = -kx*ky;
      LD[n1][n2][1][0] = -kx*ky;
      LD[n1][n2][1][1] = -ky*ky;

      LR[n1][n2][0][0] = ky*ky;
      LR[n1][n2][0][1] = -kx*ky;
      LR[n1][n2][1][0] = -kx*ky;
      LR[n1][n2][1][1] = kx*kx;
    }
  }
}