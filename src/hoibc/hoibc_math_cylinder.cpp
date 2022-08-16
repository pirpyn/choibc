#include "hoibc_math_cylinder.hpp"
#include "hoibc_math.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_bessel.hpp"
#include <cmath> // std:: pow, std::sqrt
#include <iostream>

namespace hoibc {
  namespace cylinder {
    // The tangential part of E / n x H writes in each layer where C1,C2 depends on the layer
    // E_t(r,theta,z)     = 1/(2 pi) sum_n int_R AE(r,n,kz) C1(n,kz) + BE(r,n,kz) C2(n,kz) dkz
    // n x H_t(r,theta,z) = 1/(2 pi) sum_n int_R AH(r,n,kz) C1(n,kz) + BH(r,n,kz) C2(n,kz) dkz

    // See hoibc_math_plane.cpp for more

    // Currently each cylindric computation will share the same exact impedance / reflexion.
    // As a way to speed up computation when several cylindric hoibc are computed,
    // theses will store the reflexion / impedance after the first computation.
    static big_matrix<complex> g_impedance_infinite;
    static big_matrix<complex> g_reflexion_infinite;

    matrix<complex> AE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
      matrix<complex> AE;
      const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
      AE[0][0] = -n*kz/(radius*std::pow(k3,2))*bessel1(n,k3*radius);
      AE[1][0] = bessel1(n,k3*radius);
      AE[0][1] = ci*etar*k/k3*bessel1p(n,k3*radius);
      AE[1][1] = complex(0.,0.);
      return AE;
    }

    matrix<complex> BE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
      matrix<complex> BE;
      const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
      BE[0][0] = -n*kz/(radius*std::pow(k3,2))*bessel2(n,k3*radius);
      BE[1][0] = bessel2(n,k3*radius);
      BE[0][1] = ci*etar*k/k3*bessel2p(n,k3*radius);
      BE[1][1] = complex(0.,0.);
      return BE;
    }

    matrix<complex> AH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
      matrix<complex> AH;
      const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
      AH[0][0] = complex(0.,0.);
      AH[1][0] = -ci/etar*k/k3*bessel1p(n,k3*radius);
      AH[0][1] = -bessel1(n,k3*radius);
      AH[1][1] = -n*kz/(radius*std::pow(k3,2))*bessel1(n,k3*radius);
      return AH;
    }

    matrix<complex> BH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar){
      matrix<complex> BH;
      const complex k3 = std::sqrt(std::pow(k,2) - std::pow(kz,2));
      BH[0][0] = complex(0.,0.);
      BH[1][0] = -ci/etar*k/k3*bessel2p(n,k3*radius);
      BH[0][1] = -bessel2(n,k3*radius);
      BH[1][1] = -n*kz/(radius*std::pow(k3,2))*bessel2(n,k3*radius);
      return BH;
    }

    matrix<complex> MA(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
      return AE(radius,n,kz,k,etar) - B*AH(radius,n,kz,k,etar);
    }

    matrix<complex> MB(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
      return BE(radius,n,kz,k,etar) - B*BH(radius,n,kz,k,etar);
    }

    matrix<complex> NE(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
      return AE(radius,n,kz,k,etar) + BE(radius,n,kz,k,etar)*B;
    }

    matrix<complex> NH(const real& radius,const real& n, const real& kz, const complex& k, const complex& etar, const matrix<complex>& B){
      return AH(radius,n,kz,k,etar) + BH(radius,n,kz,k,etar)*B;
    }

    big_matrix<complex> impedance_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius){
      // Was the impedance previously computed ?
      if (g_impedance_infinite.size()>0){
        return g_impedance_infinite;
      }
      big_matrix<complex> impedance_ex = big_init(vn.size(),vkz.size(),material.initial_impedance);

      if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
        std::cerr << "error: impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
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
          for (unsigned int j=0; j<vkz.size();j++) {
            const real kz = vkz[j];

            const matrix<complex>& B = impedance_ex[i][j];

            const matrix<complex> mMJ = inv(MA(radius,n,kz,k,etar,B));
            const matrix<complex> mMH = inv(MB(radius,n,kz,k,etar,B));

            const matrix<complex> mAE = AE(radius+d,n,kz,k,etar);
            const matrix<complex> mBE = BE(radius+d,n,kz,k,etar);

            const matrix<complex> mAH = AH(radius+d,n,kz,k,etar);
            const matrix<complex> mBH = BH(radius+d,n,kz,k,etar);

            // Using the overloaded / operator A / B = A*B^{-1}
            impedance_ex[i][j] = (mBE*mMH - mAE*mMJ) / (mBH*mMH - mAH*mMJ);
          }
        }
        radius = radius + d;
      }
      g_impedance_infinite = impedance_ex;
      return impedance_ex;
    }


    big_matrix<complex> reflexion_from_impedance(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> impedance, const real& outer_radius){
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
            reflexion[i][j] = - MB(radius,n,kz,k,etar,impedance[i][j]) % MA(radius,n,kz,k,etar,impedance[i][j]);
          }
        }
      }
      return reflexion;
    }

    big_matrix<complex> reflexion_infinite(const array<real>& vn, const array<real>& vkz, const real& k0, const material_t& material, const real& inner_radius){
      // Was the reflexion computed ?
      if (g_reflexion_infinite.size()>0){
        return g_reflexion_infinite;
      }
      big_matrix<complex> reflexion_ex = big_init(vn.size(),vkz.size(),complex(0.,0.));

      if (!((material.thickness.size()==material.epsr.size()) && (material.epsr.size() == material.mur.size()))){
        std::cerr << "error: impedance_infinite: size(thickness)<>size(epsr)<>size(mur)" << std::endl;
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
            - MB(radius,n,kz,k_upper,etar_upper,material.initial_impedance)
            % MA(radius,n,kz,k_upper,etar_upper,material.initial_impedance);
        }
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
          for (unsigned int j=0; j<vkz.size();j++) {
            const real kz = vkz[j];

            const matrix<complex>& B = reflexion_ex[i][j];

            const matrix<complex> mNE = inv(NE(radius,n,kz,k_lower,etar_lower,B));
            const matrix<complex> mNH = inv(NH(radius,n,kz,k_lower,etar_lower,B));

            const matrix<complex> mAE = AE(radius,n,kz,k_upper,etar_upper);
            const matrix<complex> mBE = BE(radius,n,kz,k_upper,etar_upper);

            const matrix<complex> mAH = AH(radius,n,kz,k_upper,etar_upper);
            const matrix<complex> mBH = BH(radius,n,kz,k_upper,etar_upper);

            reflexion_ex[i][j] = - ( mNE*mBE - mNH*mBH ) % ( mNE*mAE - mNH*mAH );
          }
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
        for (unsigned int j=0; j<vkz.size();j++) {
          const real kz = vkz[j];

          //if (std::abs(std::pow(k_upper,2) - std::pow(kz,2)) > 0.) {
            const matrix<complex>& B = reflexion_ex[i][j];

            const matrix<complex> mNE = inv(NE(radius,n,kz,k_lower,etar_lower,B));
            const matrix<complex> mNH = inv(NH(radius,n,kz,k_lower,etar_lower,B));

            const matrix<complex> mAE = AE(radius,n,kz,k_upper,etar_upper);
            const matrix<complex> mBE = BE(radius,n,kz,k_upper,etar_upper);

            const matrix<complex> mAH = AH(radius,n,kz,k_upper,etar_upper);
            const matrix<complex> mBH = BH(radius,n,kz,k_upper,etar_upper);

            reflexion_ex[i][j] = - ( mNE*mBE - mNH*mBH ) % ( mNE*mAE - mNH*mAH );

    //      } else {
    //        std::cerr << "warning: reflexion_infinite_cylinder: can't compute Hankel function at k3*outer_radius for k3 = 0." << std::endl;
    //        R(n1,n2,:,:) = reshape(cmplx([-1.,0.,0.,1.],[0.,0.,0.,0.],wp),[2,2])
    //      }
        }
      }
      g_reflexion_infinite = reflexion_ex;
      return reflexion_ex;
    }

    big_matrix<complex> impedance_from_reflexion(const array<real>& vn, const array<real>& vkz, const real& k0, const big_matrix<complex> reflexion, const real& outer_radius){
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

    void get_matrices_LD_LR(const real& radius, const array<real>& vn, const array<real>& vkz, big_matrix<real>& LD, big_matrix<real>& LR){
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

    void get_matrices_AB_EH(const real& radius, const array<real>& vn, const array<real>& vkz, const complex& k, const complex& etar, big_matrix<complex>& mAE ,big_matrix<complex>& mBE ,big_matrix<complex>& mAH ,big_matrix<complex>& mBH){
      mAE = big_init(vn.size(),vkz.size(),complex(0.,0.));
      mBE = big_init(vn.size(),vkz.size(),complex(0.,0.));
      mAH = big_init(vn.size(),vkz.size(),complex(0.,0.));
      mBH = big_init(vn.size(),vkz.size(),complex(0.,0.));
      for (std::size_t i = 0; i < vn.size(); i++){
        real n = vn[i];
        for (std::size_t j = 0; j < vkz.size(); j++){
          real kz = vkz[j];
          mAE[i][j] = AE(radius,n,kz,k,etar);
          mAH[i][j] = AH(radius,n,kz,k,etar);
          mBE[i][j] = BE(radius,n,kz,k,etar);
          mBH[i][j] = BH(radius,n,kz,k,etar);
        }
      }
    }

  }
}