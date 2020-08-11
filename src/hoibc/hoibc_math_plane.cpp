#include "hoibc_math_plane.hpp"
#include "hoibc_types.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_math.hpp" // matmul, inv
#include <iostream>
#include <numeric> // std::accumulate
#include <cmath> // pow, sqrt, imag, abs

using namespace hoibc;

matrix<complex> hoibc::AE(const real& kx, const real& ky, const complex& k, const real& z){
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
  AE[0][1] = AE[0][1];
  AE[1][1] = AE[0][0];
  return AE;
}

matrix<complex> hoibc::BE(const real& kx, const real& ky, const complex& k, const real& z){
  complex k3;
  if ((z>=0.)&&(std::abs(std::imag(k))<=0.)) {
    // See AE above
    k3 = std::sqrt(complex(std::pow(std::real(k),2) - std::pow(kx,2) - std::pow(ky,2),-0.));
  }
  else {
    k3 = std::sqrt(std::pow(k,2) - std::pow(kx,2) - std::pow(ky,2));
  }
  matrix<complex> BE;
  BE[0][0] = std::exp(-ci*k3*z)*ci*k3;
  BE[1][0] = 0.;
  BE[0][1] = BE[1][0];
  BE[1][1] = BE[0][0];
  return BE;
}

matrix<complex> hoibc::AH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z){
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
  AH[1][0] = std::exp(ci*k3*z)*ci/k/etar*(-kx*ky);
  AH[0][1] = AH[1][0];
  AH[1][1] = std::exp(ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(kx,2));
  return AH;
}

matrix<complex> hoibc::BH(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z){
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
  BH[1][0] = std::exp(-ci*k3*z)*ci/k/etar*(-kx*ky);
  BH[0][1] = BH[1][0];
  BH[1][1] = std::exp(-ci*k3*z)*ci/k/etar*(std::pow(k,2)-std::pow(kx,2));
  return BH;
}

matrix<complex> hoibc::MA(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp){
  return AE(kx,ky,k,z) - matmul(imp,AH(kx,ky,k,etar,z));
}

matrix<complex> hoibc::MB(const real& kx, const real& ky, const complex& k, const complex& etar, const real& z, const matrix<complex>& imp){
  return BE(kx,ky,k,z) - matmul(imp,BH(kx,ky,k,etar,z));
}

big_matrix<complex> hoibc::impedance_infinite_plane(const std::vector<real> &vkx, const  std::vector<real> &vky, const real& k0, const material_t& material){

  big_matrix<complex> impedance_ex = big_init(vkx.size(), vky.size(), material.initial_impedance);

  const std::vector<real>& thickness = material.thickness;

  real h = - std::accumulate(thickness.begin(),thickness.end(),static_cast<real>(0));

  for (unsigned int i=0;i<thickness.size();i++) {
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

    const complex k = k0*etar;
    for (unsigned int i=0; i<vkx.size();i++) {
      const real kx = vkx[i];
      for (unsigned int j=0; j<vky.size();j++) {
        const real ky = vky[j];

        const matrix<complex>& B = impedance_ex[i][j];

        const matrix<complex> mBE = BE(kx,ky,k,h+d);
        const matrix<complex> mAE = AE(kx,ky,k,h+d);
        const matrix<complex> mBH = BH(kx,ky,k,etar,h+d);
        const matrix<complex> mAH = AH(kx,ky,k,etar,h+d);
        const matrix<complex> mMB = inv<complex>(MB(kx,ky,k,etar,h,B));
        const matrix<complex> mMA = inv<complex>(MA(kx,ky,k,etar,h,B));

        impedance_ex[i][j] = matmul<complex>(
          matmul<complex>(mBE,mMB) - matmul<complex>(mAE,mMA),
          inv<complex>(matmul<complex>(mBH,mMB) - matmul<complex>(mAH,inv<complex>(mMA)))
          );
      }
    }
    h+=d;
  }
  return impedance_ex;
}

void hoibc::reflexion_infinite_plane(){
  std::cout << "hoibc_math_plane::reflexion_infinite_plane: i do nothing" << std::endl;
}
