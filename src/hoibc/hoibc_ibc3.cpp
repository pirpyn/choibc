#include "hoibc_ibc3.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"

#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include "lapacke.h"

using namespace hoibc;

void hoibc::hoibc_ibc3::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  
  big_matrix<real> I;
  const std::size_t n1 = f1.size();
  const std::size_t n2 = f2.size();

  get_matrices_I(n1,n2,I=I);
  
  // Theses are used when this->type = 1
  big_matrix<complex> AE, AH, BE, BH;
  big_matrix<complex> ME, MH;
  big_matrix<complex> LDME, LDMH;
  big_matrix<complex> LRME, LRMH;

  // Theses are used when this->type = 2
  big_matrix<complex> LDZ;
  big_matrix<complex> LRZ;
  
  big_matrix<real> LD, LR;
  this->get_matrices_LD_LR(k0,f1,f2,LD,LR);

  lapack_int m = f1.size()*f2.size()*4;
  lapack_int n = 5;
  // The array to pass to LAPACK that contains the matrix in the 
  // least square resolution and its rhs.
  complex **a = new complex*[m];
  a[0] = new complex[m * n];
  for (int i = 1; i < m; i++)
      a[i] = a[i-1] + n;
  complex* b = new complex[m];

  switch (this->mode){
  case mode_t::R: // if reflexion
    this->get_matrices_AB_EH(k0,f1,f2,AE,BE,AH,BH);
    ME = AE + BE * gex;
    MH = AH + BH * gex;
    LDME = LD * ME;
    LDMH = LD * MH;
    LRME = LR * ME;
    LRMH = LR * MH;
    for (std::size_t j=0;j<n2;j++){
      for (std::size_t i=0;i<n1;i++){
        a[(j*n1+i)*4+0][0] = MH[i][j][0][0];
        a[(j*n1+i)*4+1][0] = MH[i][j][1][0];
        a[(j*n1+i)*4+2][0] = MH[i][j][0][1];
        a[(j*n1+i)*4+3][0] = MH[i][j][1][1];

        a[(j*n1+i)*4+0][1] = LDMH[i][j][0][0];
        a[(j*n1+i)*4+1][1] = LDMH[i][j][1][0];
        a[(j*n1+i)*4+2][1] = LDMH[i][j][0][1];
        a[(j*n1+i)*4+3][1] = LDMH[i][j][1][1];

        a[(j*n1+i)*4+0][2] = -(LRMH[i][j][0][0]);
        a[(j*n1+i)*4+1][2] = -(LRMH[i][j][1][0]);
        a[(j*n1+i)*4+2][2] = -(LRMH[i][j][0][1]);
        a[(j*n1+i)*4+3][2] = -(LRMH[i][j][1][1]);

        a[(j*n1+i)*4+0][3] = -(LDME[i][j][0][0]);
        a[(j*n1+i)*4+1][3] = -(LDME[i][j][1][0]);
        a[(j*n1+i)*4+2][3] = -(LDME[i][j][0][1]);
        a[(j*n1+i)*4+3][3] = -(LDME[i][j][1][1]);

        a[(j*n1+i)*4+0][4] = LRME[i][j][0][0];
        a[(j*n1+i)*4+1][4] = LRME[i][j][1][0];
        a[(j*n1+i)*4+2][4] = LRME[i][j][0][1];
        a[(j*n1+i)*4+3][4] = LRME[i][j][1][1];

        b[(j*n1+i)*4+0] = ME[i][j][0][0];
        b[(j*n1+i)*4+1] = ME[i][j][1][0];
        b[(j*n1+i)*4+2] = ME[i][j][0][1];
        b[(j*n1+i)*4+3] = ME[i][j][1][1];
      }
    }
    break;
  case mode_t::Z: // if impedance
    LDZ = LD * gex;
    LRZ = LR * gex;

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

        a[(j*n1+i)*4+0][2] = -(LR[i][j][0][0]);
        a[(j*n1+i)*4+1][2] = -(LR[i][j][1][0]);
        a[(j*n1+i)*4+2][2] = -(LR[i][j][0][1]);
        a[(j*n1+i)*4+3][2] = -(LR[i][j][1][1]);

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
    break;
  default:
    std::cerr << "hoibc_ibc3::get_coeff_no_suc: mode = " << mode_to_int(this->mode) << " unknown" << std::endl;
    exit(1);
    break;
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
}

big_matrix<complex> hoibc::hoibc_ibc3::get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2){
  big_matrix<real> I, LD, LR;

  get_matrices_I(f1.size(),f2.size(),I=I);

  this->get_matrices_LD_LR(k0,f1,f2,LD,LR);

  big_matrix<complex> impedance_ap = 
      ( complex(1.,0.)*I + this->coeff.b1*LD - this->coeff.b2*LR ) 
    % ( this->coeff.a0*I + this->coeff.a1*LD - this->coeff.a2*LR );
  return impedance_ap;
}

void hoibc::hoibc_ibc3::array_to_coeff(const alglib::real_1d_array& x){
  this->coeff.a0 = complex(x[0],x[1]);
  this->coeff.a1 = complex(x[2],x[3]);
  this->coeff.a2 = complex(x[4],x[5]);
  this->coeff.b1 = complex(x[6],x[7]);
  this->coeff.b2 = complex(x[8],x[9]);
}

void hoibc::hoibc_ibc3::coeff_to_array(alglib::real_1d_array& x){
  x.setlength(10);
  x[0] = std::real(this->coeff.a0);
  x[1] = std::imag(this->coeff.a0);
  x[2] = std::real(this->coeff.a1);
  x[3] = std::imag(this->coeff.a1);
  x[4] = std::real(this->coeff.a2);
  x[5] = std::imag(this->coeff.a2);
  x[6] = std::real(this->coeff.b1);
  x[7] = std::imag(this->coeff.b1);
  x[8] = std::real(this->coeff.b2);
  x[9] = std::imag(this->coeff.b2);
}

void hoibc::hoibc_ibc3::get_suc(std::vector<real>& cle, std::vector<real>& ceq, std::vector<real>& cne, std::vector<std::string>& sle, std::vector<std::string>& seq, std::vector<std::string>& sne){
    cle.clear();
    ceq.clear();
    cne.clear();
    sle.clear();
    seq.clear();
    sne.clear();

    coeff_t &c = this->coeff;
    const complex z = std::pow(std::abs(c.a1),2)*std::pow(std::abs(c.a2),2) - c.b1*c.a0*conj(c.a1)*std::pow(std::abs(c.a2),2) - c.b2*c.a0*conj(c.a2)*std::pow(std::abs(c.a1),2);

    cle.resize(9);
    cle[0] = -std::real(conj(c.a0)*z);
    cle[1] =  std::real(conj(c.a1)*z);
    cle[2] =  std::real(conj(c.a2)*z);
    cle[3] = -std::real(c.b1*conj(c.a1));
    cle[4] = -std::real(c.b2*conj(c.a2));
    cle[5] =  std::real(c.b1*conj(c.a2*c.a1)*c.a0);
    cle[6] =  std::real(c.b2*conj(c.a1*c.a2)*c.a0);
    cle[7] =  std::real(c.a1);
    cle[8] =  std::real(c.a2);
    sle.resize(9);
    sle[0] = "-Re(a0*.z); z = |a1|^2.|a2|^2 - b1.a0.a1*.|a2|^2 - b2.a0.a2*.|a1|^2";
    sle[1] = " Re(a1*.z); z = |a1|^2.|a2|^2 - b1.a0.a1*.|a2|^2 - b2.a0.a2*.|a1|^2";
    sle[2] = " Re(a2*.z); z = |a1|^2.|a2|^2 - b1.a0.a1*.|a2|^2 - b2.a0.a2*.|a1|^2";
    sle[3] = "-Re(b1.a1*)";
    sle[4] = "-Re(b2.a2*)";
    sle[5] = " Re(b1.(a2.a1)*.a0)";
    sle[6] = " Re(b2.(a2.a1)*.a0)";
    sle[7] = " Re(a1)";
    sle[8] = " Re(a2)";
}

void hoibc::hoibc_ibc3::disp_coeff(std::ostream& out){
  if (this->normalised){
    out << "# Z = (I + b1*LD/k0^2 - b2*LR/k0^2)^{-1} (a0*I + a1*LD/k0^2 - a2*LR/k0^2)" << std::endl;
  } else {
    out << "# Z = (I + b1*LD - b2*LR)^{-1} (a0*I + a1*LD - a2*LR)" << std::endl;
  }
  
  print_complex(this->coeff.a0,"  a0",out);
  print_complex(this->coeff.a1,"  a1",out);
  print_complex(this->coeff.a2,"  a2",out);
  print_complex(this->coeff.b1,"  b1",out);
  print_complex(this->coeff.b2,"  b2",out);
}

void hoibc::hoibc_ibc3::get_matrices_LD_LR(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2, big_matrix<real>& LD, big_matrix<real>& LR){
    switch (this->type){
      case type_t::P:
        plane::get_matrices_LD_LR(f1,f2,LD=LD,LR=LR);
        break;
      default:
        NOTFINISHED("hoibc::hoibc_ibc3::get_matrices_LD_LR")
      break;
    }
    if (this->normalised){
      LD = LD / (k0*k0);
      LR = LR / (k0*k0);
    }
}

void hoibc::hoibc_ibc3::get_matrices_AB_EH(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2, big_matrix<complex>& AE, big_matrix<complex>& BE, big_matrix<complex>& AH, big_matrix<complex>& BH){

  // We're in the vaccum
  const complex k = k0;
  const complex etar = 1.;
  const real height = 0.;

  switch (this->type){
  case type_t::P:
    plane::get_matrices_AB_EH(f1,f2,k,etar,height,AE,BE,AH,BH);
    break;
  // case ('C'):
  //   cylinder::get_matrices_AB_EH(this->outer_radius,f1,f2,k,etar,AE,BE,AH,BH);
  //   break;
  // case ('S'):
  //   sphere::get_matrices_AB_EH(this->outer_radius,f2,k,etar,AE,BE,AH,BH);
  //   break;
  default:
    NOTFINISHED("hoibc::hoibc_ibc3::get_matrices_AB_EH")
    break;
  }
}