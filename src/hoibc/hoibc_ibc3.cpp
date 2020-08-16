#include "hoibc_ibc3.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"

#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include "lapacke.h"

using namespace hoibc;

void hoibc::hoibc_ibc3::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  
  big_matrix<real> I, LD, LR;
  const std::size_t n1 = f1.size();
  const std::size_t n2 = f2.size();
  complex** a;
  complex* b;
  lapack_int m, n;
  get_matrices_I(n1,n2,I=I);
  big_matrix<complex> LDZ;
  big_matrix<complex> LRZ;
  switch (this->mode){
  // case 1: // if reflexion
  //   NOTFINISHED("hoibc::hoibc_ibc3_get_coeff_no_suc::mode=1")
  //  break;
  case 2: // if impedance
    switch (this->type){
      case 'P':
        plane::get_matrices_LD_LR(f1,f2,LD=LD,LR=LR);
        break;
      default:
        NOTFINISHED("hoibc::ibc3::get_coeff_no_suc::mode=2")
      break;
    }
    if (this->normalised){
      LD = LD / (k0*k0);
      LR = LR / (k0*k0);
    }
    LDZ = LD * gex;
    LRZ = LR * gex;
    m = f1.size()*f2.size()*4;
    n = 5;

    // Allocating a contiguous memory for the "2D" array to benefits from LAPACK efficiency.
    a = new complex*[m];
    a[0] = new complex[m * n];
    for (int i = 1; i < m; i++)
        a[i] = a[i-1] + n;
    b = new complex[m];

    for (std::size_t j=0;j<n2;j++){
      for (std::size_t i=0;i<n1;i++){
        a[(j*n1+i)*4+0][0] = 1.;
        a[(j*n1+i)*4+1][0] = 0.;
        a[(j*n1+i)*4+2][0] = 0.;
        a[(j*n1+i)*4+3][0] = 1.;

        a[(j*n1+i)*4+0][1] = LD[i][j][0][0];
        a[(j*n1+i)*4+1][1] = LD[i][j][1][0];
        a[(j*n1+i)*4+2][1] = LD[i][j][0][1];
        a[(j*n1+i)*4+3][1] = LD[i][j][1][1];

        a[(j*n1+i)*4+0][2] = -LR[i][j][0][0];
        a[(j*n1+i)*4+1][2] = -LR[i][j][1][0];
        a[(j*n1+i)*4+2][2] = -LR[i][j][0][1];
        a[(j*n1+i)*4+3][2] = -LR[i][j][1][1];

        a[(j*n1+i)*4+0][3] = -(LDZ[i][j][0][0]);
        a[(j*n1+i)*4+1][3] = -(LDZ[i][j][1][0]);
        a[(j*n1+i)*4+2][3] = -(LDZ[i][j][0][1]);
        a[(j*n1+i)*4+3][3] = -(LDZ[i][j][1][1]);

        a[(j*n1+i)*4+0][4] = LRZ[i][j][0][0];
        a[(j*n1+i)*4+1][4] = LRZ[i][j][1][0];
        a[(j*n1+i)*4+2][4] = LRZ[i][j][0][1];
        a[(j*n1+i)*4+3][4] = LRZ[i][j][1][1];

        b[(j*n1+i)*4+0] = gex[i][j][0][0];
        b[(j*n1+i)*4+1] = gex[i][j][1][0];
        b[(j*n1+i)*4+2] = gex[i][j][0][1];
        b[(j*n1+i)*4+3] = gex[i][j][1][1];
      }
    }

    LAPACKE_zgels(LAPACK_ROW_MAJOR,'N', m, n, 1, *a, n, b, 1 );

    this->coeff.a0 = b[0];
    this->coeff.a1 = b[1];
    this->coeff.a2 = b[2];
    this->coeff.b1 = b[3];
    this->coeff.b2 = b[4];
    delete [] a[0];
    delete [] a;
    delete [] b;
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

  big_matrix<complex> impedance_ap = ( complex(1.,0.)*I + (this->coeff.b1*LD) - (this->coeff.b2*LR) ) % ( this->coeff.a0*I + this->coeff.a1*LD - this->coeff.a2*LR );
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
  if (this->normalised){
    out << "# Z = (I + b1*LD/k0^2 - b2*LR^2)^{-1}*(a0*I + a1*LD^2 - a2*LR^2)" << std::endl;
  } else {
    out << "# Z = (I + b1*LD - b2*LR)^{-1}*(a0*I + a1*LD - a2*LR)" << std::endl;
  }
  
  print_complex(this->coeff.a0,"    a0",out);
  print_complex(this->coeff.a1,"    a1",out);
  print_complex(this->coeff.a2,"    a2",out);
  print_complex(this->coeff.b1,"    b1",out);
  print_complex(this->coeff.b2,"    b2",out);
}