#include "hoibc_bessel.hpp"
#include "../bessel/zbessel.hh"
#include "hoibc_constants.hpp"
#include <cassert>
#include <string>
#include <iostream>

namespace hoibc {
  static std::string error_msg(const int& ierr, const real& nu, const complex& z){
    switch (ierr){
    case 1:
      return  "Input Error        - NO COMPUTATION";
    case 2:
      return  "Overflow           - NO COMPUTATION (Im(Z) too large with KODE=1)"
              "                     Im(Z) = "+std::to_string(z.imag());
    case 3:
      return  "Precision warning  - COMPUTATION COMPLETED (Result has half precision or less because abs(z) or FNU+N-1 is large)"
              "                     Abs(Z) = "+std::to_string(std::abs(z))+", nu = "+std::to_string(nu);
    case 4:
      return  "Precision error    - NO COMPUTATION (Result has no precision because abs(z) or FNU+N-1 is too large)"
              "                     Abs(Z) = "+std::to_string(std::abs(z))+", nu = "+std::to_string(nu);
    case 5:
      return  "Algortithmic error - NO COMPUTATION (Termination condition not met)";
    }
    return    "Normal return      - COMPUTATION COMPLETED";
  }

  complex bessel1(const real& nu, const complex& z){
    real yr, yi;
    int nz;
    const int ierr = zbessel::zbesj(z.real(),z.imag(),nu,1,1,&yr,&yi,&nz);
    switch(ierr){
    case 1:
    case 5:
      std::cerr << "error: zbesj(nu,z), nu = " << nu << ",z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      exit(ierr);
      break;
    case 3:
      std::cerr << "warning: zbesj(nu,z), nu = " << nu << ", z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      break;
    case 2:
    case 4:
      std::cerr << "critical warning: with zbesj(nu,z), nu = " << nu << ", z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      break;
    }
    return complex(yr,yi);
  }

  complex bessel2(const real& nu, const complex& z){
    real yr, yi;
    int nz;
    const int ierr = zbessel::zbesh(z.real(),z.imag(),nu,1,2,1,&yr,&yi,&nz);
    switch(ierr){
    case 1:
    case 5:
      std::cerr << "error: zbesj(nu,z), nu = " << nu << ",z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      exit(ierr);
      break;
    case 3:
      std::cerr << "warning: zbesj(nu,z), nu = " << nu << ", z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      break;
    case 2:
    case 4:
      std::cerr << "critical warning: with zbesj(nu,z), nu = " << nu << ", z = " << z << std::endl;
      std::cerr << error_msg(ierr,nu,z) << std::endl;
      break;
    }
    return complex(yr,yi);
  }

  complex bessel1p(const real& nu, const complex& z){
    return - bessel1(nu+1,z) + nu/z*bessel1(nu,z);
  }

  complex bessel2p(const real& nu, const complex& z){
    return - bessel2(nu+1,z) + nu/z*bessel2(nu,z);
  }

  complex sbessel1(const real& nu, const complex& z){
    return std::sqrt(0.5*pi/z)*bessel1(nu+0.5,z);
  }

  complex sbessel2(const real& nu, const complex& z){
    return std::sqrt(0.5*pi/z)*bessel2(nu+0.5,z);
  }

  complex sbessel1p(const real& nu, const complex& z){
    return std::sqrt(0.5*pi/z)*bessel1p(nu+0.5,z) - 0.5*std::sqrt(0.5*pi/std::pow(z,3))*bessel1(nu+0.5,z);
  }

  complex sbessel2p(const real& nu, const complex& z){
    return std::sqrt(0.5*pi/z)*bessel2p(nu+0.5,z) - 0.5*std::sqrt(0.5*pi/std::pow(z,3))*bessel2(nu+0.5,z);
  }

}