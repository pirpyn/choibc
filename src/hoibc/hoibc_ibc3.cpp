#include "hoibc_ibc3.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"

#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include "lapacke.h"

using namespace hoibc;

void hoibc::hoibc_ibc3::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  // We're in the vaccum
  const complex k = k0;
  const complex etar = 1.;
  const real z = 0.;

  big_matrix<complex> mAE,mBE,mAH,mBH;
  big_matrix<real> I;

  switch (this->mode){
  case 1: // if reflexion
    NOTFINISHED("hoibc::hoibc_ibc3_get_coeff_no_suc::mode=1")
  case 2: // if impedance
      const std::size_t n1 = f1.size();
      const std::size_t n2 = f2.size();
      get_matrices_I(n1,n2,I=I);

      lapack_int m = f1.size()*f2.size()*4;// std::vector<std::array<complex,5>> A;
      lapack_int n = 5;
      complex** a = new complex[m][n];
      complex** b = new complex[m];

      for (std::size_t j=0;j<n2;j++){
        for (std::size_t i=0;i<n1;i++){
          a[(j*n1+i)*4+0][0] = 1.;
          a[(j*n1+i)*4+1][0] = 0.;
          a[(j*n1+i)*4+2][0] = 0;
          a[(j*n1+i)*4+3][0] = 1.;
        }
      }
    //   A(:,1) = reshape(I,[iM])
    //   A(:,2) = reshape(LD,[iM])
    //   A(:,3) = reshape(-LR,[iM])
    //   A(:,4) = reshape(-big_matmul(LD,exact),[iM])
    //   A(:,5) = reshape(big_matmul(LR,exact),[iM])

    //   b(1:iM) = reshape(exact,[iM])
    // end select

    // LAPACKE_zgels(LAPACK_COL_MAJOR,'N', m, n, 1, *a, m, *b, m );

    // this->coeff.a0 = b[0];
    // this->coeff.a1 = b[1];
    // this->coeff.a2 = b[2];
    // this->coeff.b1 = b[3];
    // this->coeff.b2 = b[4];
    delete a;
    delete b;
    break;
  default:
    std::cerr << "hoibc_ibc3::get_coeff_no_suc: mode = " << this->mode << " unknown" << std::endl;
    break;
  }
  NOTFINISHED("hoibc_ibc3::get_coeff_no_suc")
}

big_matrix<complex> hoibc::hoibc_ibc3::get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2){
  big_matrix<real> I, LD, LR;

  get_matrices_I(f1.size(),f2.size(),I=I);

  switch (this->type){
    case 'P':
      plane::get_matrices_LD_LR(f1,f2,LD=LD,LR=LR);
      break;
    default:
      NOTFINISHED("hoibc_ibc3::get_impedance")
      break;
  }

  if (this->normalised){
    const real scaling = k0*k0;
    LD = LD / scaling;
    LR = LR / scaling;
  }

  big_matrix<complex> impedance_ap = ( complex(1.,0.)*I + (this->coeff.b1*LD) + (this->coeff.b2*LR) ) % ( this->coeff.a0*I + this->coeff.a1*LD + this->coeff.a2*LR );
  return impedance_ap;
}

void hoibc::hoibc_ibc3::array_to_coeff(const std::vector<real>& x){
  this->coeff.a0 = std::complex<real>(x[0],x[1]);
  this->coeff.a1 = std::complex<real>(x[2],x[3]);
  this->coeff.a2 = std::complex<real>(x[4],x[5]);
  this->coeff.b1 = std::complex<real>(x[6],x[7]);
  this->coeff.b2 = std::complex<real>(x[8],x[9]);
}

void hoibc::hoibc_ibc3::coeff_to_array(std::vector<real>& x){
  x.clear();
  x.push_back(std::real(this->coeff.a0));
  x.push_back(std::imag(this->coeff.a0));
  x.push_back(std::real(this->coeff.a1));
  x.push_back(std::imag(this->coeff.a1));
  x.push_back(std::real(this->coeff.a2));
  x.push_back(std::imag(this->coeff.a2));
  x.push_back(std::real(this->coeff.b1));
  x.push_back(std::imag(this->coeff.b1));
  x.push_back(std::real(this->coeff.b2));
  x.push_back(std::imag(this->coeff.b2));
}

void hoibc::hoibc_ibc3::get_suc(std::vector<real>& cle, std::vector<real>& ceq, std::vector<real>& cne, std::vector<std::string>& sle, std::vector<std::string>& seq, std::vector<std::string>& sne){
  cle.clear();
  ceq.clear();
  cne.clear();
  sle.clear();
  seq.clear();
  sne.clear();
}

void hoibc::hoibc_ibc3::disp_coeff(std::ostream& out){
  out << "# Z = (I + b1*LD - b2*LR)^{-1}*(a0*I + a1*LD - a2*LR)" << std::endl;
  print_complex(this->coeff.a0,"    a0",out);
  print_complex(this->coeff.a1,"    a1",out);
  print_complex(this->coeff.a2,"    a2",out);
  print_complex(this->coeff.b1,"    b1",out);
  print_complex(this->coeff.b2,"    b2",out);
}